module SolveLAP
##### PLEASE READ THE NOTICE ON THE BOTTOM, for the original license and authors
using SparseArrays
export solve_lap, solve_stiff_lap
include("StiffWrapper.jl")
# Zero markers used in hungarian algorithm
# 0 => NON   => Non-zero
# 1 => Z     => ordinary zero
# 2 => STAR  => starred zero
# 3 => PRIME => primed zero
const NON = Int8(0)
const Z = Int8(1)
const STAR = Int8(2)
const PRIME = Int8(3)

"""
solve_lap(cost_matrix) -> (assignment, cost)

# Example
```julia
julia> using SolveLAP

julia> M=rand(1:100,3,3)
3×3 Array{Int64,2}:
 66  83  92
 17  92  65
 30  84  18

julia> solve_lap(M)
([2, 1, 3], 118)
```
"""
function solve_lap(costMat::AbstractMatrix)
    rowNum, colNum = size(costMat)
    # currently, the function `hungarian` automatically transposes `cost matrix` when there are more workers than jobs.
    costMatrix = rowNum <= colNum ? costMat : transpose(costMat)
    matching = build_matching(costMatrix)
    assignment = zeros(Int, rowNum)
    rows = rowvals(matching)
    for c = 1:size(matching,2), i in nzrange(matching, c)
        r = rows[i]
        if matching[r,c] == STAR
            if rowNum ≤ colNum
                assignment[r] = c
            else
                assignment[c] = r
            end
        end
    end
    # calculate minimum cost
    cost = sum(costMat[i...] for i in zip(1:rowNum, assignment) if i[2] != 0)
    return assignment,cost
end

"""
solve_stiff_lap(cost_matrix) -> assignment in cartesian coords

# Example
```julia
julia> using SolveLAP

julia> M=rand([1,2,Inf],5,5)
5×5 Array{Float64,2}:
 Inf    Inf      1.0    1.0    1.0
 Inf    Inf    Inf    Inf    Inf
   2.0    1.0  Inf      1.0    2.0
 Inf      2.0    1.0    2.0  Inf
   2.0    2.0    2.0  Inf      1.0

julia> solve_stiff_lap(M)
6-element Array{Tuple{Int64,Int64},1}:
 (1, 4)        ← row 1 is assigned to column 4
 (2, -1)       ← -1 for missing assignment, row 2 is assigned to nothing
 (3, 2)
 (4, 3)
 (5, 5)
 (-1, 1)       ← column 1 is assigned to nothing
```
"""
function solve_stiff_lap(C::AbstractMatrix{T}, penalty=1.05) where T <: Real
    m,n=size(C)
    assignment,cost=solve_lap(StiffWrapper(C,penalty))
    c=[(i,assignment[i]) for i in 1:(m+n)]
    for k in eachindex(c)
        i,j=c[k]
        i<=m||(i=-1)
        j<=n||(j=-1)
        c[k]=(i,j)
    end
    [a for a in c if a[1]+a[2]>-2]
end

function build_matching(costMat::AbstractMatrix{T}) where T <: Real
    rowNum, colNum = size(costMat)
    colNum >= rowNum || throw(ArgumentError("Non-square matrix should have more columns than rows."))

    # preliminaries:
    # "no lines are covered;"
    rowCovered = falses(rowNum)
    colCovered = falses(colNum)


    # for tracking indices
    rowCoveredIdx = Int[]
    colCoveredIdx = Int[]
    rowUncoveredIdx = Int[]
    colUncoveredIdx = Int[]
    sizehint!(rowCoveredIdx, rowNum)
    sizehint!(colCoveredIdx, colNum)
    sizehint!(rowUncoveredIdx, rowNum)
    sizehint!(colUncoveredIdx, colNum)

    # "no zeros are starred or primed."
    # use a sparse matrix Zs to store these three kinds of zeros:
    # 0 => NON   => Non-zero
    # 1 => Z     => ordinary zero
    # 2 => STAR  => starred zero
    # 3 => PRIME => primed zero
    Zs = spzeros(Int8, rowNum, colNum)

    # for tracking changes per row/col of A
    Δrow = zeros(T,rowNum)
    Δcol = zeros(T,colNum)

    # "subtract" minimum from each row
    @inbounds for i in 1:rowNum
        mw=costMat[i,1]
        @inbounds for j in 2:colNum
            cost=costMat[i,j]
            if cost<mw
                mw=cost
            end
        end
        Δrow[i] = -mw
    end
    # "subtract" minimum from each column
    if rowNum == colNum
        @inbounds for j in 1:colNum
            mw=costMat[1,j]+Δrow[1]
            @inbounds for i in 2:rowNum
                cost=costMat[i,j]+Δrow[i]
                if cost<mw
                    mw=cost
                end
            end
            Δcol[j] = -mw
        end
    end
    for i in Δcol
        if isinf(i)
            error("A column contains only Inf, unsolvable problem via solve_lap, try solve_stiff_lap")
        end
    end
    for i in Δrow
        if isinf(i)
            error("A row contains only Inf, unsolvable problem via solve_lap, try solve_stiff_lap")
        end
    end

    # for tracking those starred zero
    rowSTAR = falses(rowNum)
    colSTAR = falses(colNum)
    # since looping through a row in a SparseMatrixCSC is costy(not cache-friendly),
    # a row to column index mapping of starred zeros is also tracked here.
    row2colSTAR = Dict{Int,Int}()
    for ii in CartesianIndices(costMat)
        r,c=Tuple(ii)
        cost=costMat[r,c]+Δrow[r]+Δcol[c]
        # "consider a zero Z of the matrix;"
        if iszero(cost)
            Zs[r,c] = Z
            # "if there is no starred zero in its row and none in its column, star Z.
            #  repeat, considering each zero in the matrix in turn;"
            if !colSTAR[c] && !rowSTAR[r]
                Zs[r,c] = STAR
                rowSTAR[r] = true
                colSTAR[c] = true
                row2colSTAR[r] = c
                # "then cover every column containing a starred zero."
                colCovered[c] = true
            end
        end
    end

    # preliminaries done, go to step 1
    stepNum = 1

    # if the assignment is already found, exit
    stepNum = exit_criteria(colCovered, size(Zs))

    # pre-allocation
    sequence = Tuple{Int,Int}[]       # used in step 2
    minLocations = Tuple{Int,Int}[]   # used in step 3

    # these three steps are parallel with those in the paper:
    # J. Munkres, "Algorithms for the Assignment and Transportation Problems",
    # Journal of the Society for Industrial and Applied Mathematics, 5(1):32–38, 1957 March.
    while stepNum != 0
        if stepNum == 1
            stepNum = step1!(Zs, rowCovered, colCovered, rowSTAR, row2colSTAR)
        elseif stepNum == 2
            empty!(sequence)
            stepNum = step2!(Zs, sequence, rowCovered, colCovered, rowSTAR, row2colSTAR)
        elseif stepNum == 3
            empty!(rowCoveredIdx)
            empty!(colCoveredIdx)
            empty!(rowUncoveredIdx)
            empty!(colUncoveredIdx)
            empty!(minLocations)
            stepNum = step3!(costMat, Zs, minLocations, rowCovered, colCovered, rowSTAR, row2colSTAR, Δrow, Δcol, rowCoveredIdx, colCoveredIdx, rowUncoveredIdx, colUncoveredIdx)
        end
    end

    return Zs
end



"""
    exit_criteria(colCovered, ZsDims) -> stepNum
We adjust Munkres's algorithm in order to deal with rectangular matrices,
so only K columns are counted here, where K = min(size(Zs))
"""
function exit_criteria(colCovered, ZsDims)
    count = 0
    @inbounds for i in eachindex(colCovered)
        colCovered[i] && (count += 1;)
    end
    if count == ZsDims[1]
        # algorithm exits
        return 0
    else
        # "otherwise, return to step 1."
        return 1
    end
end


"""
Step 1 of the original Munkres' Assignment Algorithm
"""
function step1!(Zs, rowCovered, colCovered, rowSTAR, row2colSTAR)
    colLen = size(Zs,2)
    rows = rowvals(Zs)
    # step 1:
    zeroCoveredNum = 0
    # "repeat until all zeros are covered."
    while zeroCoveredNum < nnz(Zs)
        zeroCoveredNum = 0
        for c = 1:colLen, i in nzrange(Zs, c)
            r = rows[i]
            # "choose a non-covered zero and prime it"
            if colCovered[c] == false && rowCovered[r] == false
                Zs[r,c] = PRIME
                # "consider the row containing it."
                # "if there is a starred zero Z in this row"
                if rowSTAR[r]
                    # "cover this row and uncover the column of Z"
                    rowCovered[r] = true
                    colCovered[row2colSTAR[r]] = false
                else
                    # "if there is no starred zero in this row,
                    #  go at once to step 2."
                    return 2
                end
            else
                # otherwise, this zero is covered
                zeroCoveredNum += 1
            end
        end
    end
    # "go to step 3."
    return 3
end

"""
Step 2 of the original Munkres' Assignment Algorithm
"""
function step2!(Zs, sequence, rowCovered, colCovered, rowSTAR, row2colSTAR)
    ZsDims = size(Zs)
    rows = rowvals(Zs)
    # step 2:
    flag = false
    # "there is a sequence of alternating primed and starred zeros, constructed
    #  as follows:"
    # "let Z₀ denote the uncovered 0′.[there is only one.]"
    for c = 1:ZsDims[2], i in nzrange(Zs, c)
        r = rows[i]
        # note that Z₀ is an **uncovered** 0′
        if Zs[r,c] == PRIME && colCovered[c] == false && rowCovered[r] == false
            # push Z₀(uncovered 0′) into the sequence
            push!(sequence, (r, c))
            # "let Z₁ denote the 0* in the Z₀'s column.(if any)"
            # find 0* in Z₀'s column
            for j in nzrange(Zs, c)
                Z₁r = rows[j]
                if Zs[Z₁r, c] == STAR
                    # push Z₁(0*) into the sequence
                    push!(sequence, (Z₁r, c))
                    # set sequence continue flag => true
                    flag = true
                    break
                end
            end
            break
        end
    end
    # "let Z₂ denote the 0′ in Z₁'s row(there will always be one)."
    # "let Z₃ denote the 0* in the Z₂'s column."
    # "continue until the sequence stops at a 0′ Z₂ⱼ, which has no 0* in its column."
    while flag
        flag = false
        r = sequence[end][1]
        # find Z₂ in Z₃'s row (always exits)
        for c = 1:ZsDims[2]
            if Zs[r,c] == PRIME
                # push Z₂ into the sequence
                push!(sequence, (r, c))
                # find 0* in Z₂'s column
                for j in nzrange(Zs, c)
                    Z₃r = rows[j]
                    if Zs[Z₃r, c] == STAR
                        # push Z₃(0*) into the sequence
                        push!(sequence, (Z₃r, c))
                        # set sequence continue flag => true
                        flag = true
                        break
                    end
                end
                break
            end
        end
    end

    # "unstar each starred zero of the sequence;"
    for i in 2:2:length(sequence)-1
        Zs[sequence[i]...] = Z
    end

    # clean up
    fill!(rowSTAR, false)
    empty!(row2colSTAR)

    # "and star each primed zero of the sequence."
    for i in 1:2:length(sequence)
        Zs[sequence[i]...] = STAR
    end

    for c = 1:ZsDims[2], i in nzrange(Zs, c)
        r = rows[i]
        # "erase all primes;"
        if Zs[r,c] == PRIME
            Zs[r,c] = Z
        end
        # "and cover every column containing a 0*."
        if Zs[r,c] == STAR
            colCovered[c] = true
            rowSTAR[r] = true
            row2colSTAR[r] = c
        end
    end

    # "uncover every row"
    fill!(rowCovered, false)

    # "if all columns are covered, the starred zeros form the desired independent set."
    return exit_criteria(colCovered, ZsDims)
end

"""
Step 3 of the original Munkres' Assignment Algorithm
"""
function step3!(costMat::AbstractMatrix{T}, Zs, minLocations, rowCovered, colCovered, rowSTAR, row2colSTAR,
                Δrow, Δcol, rowCoveredIdx, colCoveredIdx, rowUncoveredIdx, colUncoveredIdx) where T <: Real

    @inbounds for i in eachindex(rowCovered)
        rowCovered[i] ? push!(rowCoveredIdx, i) : push!(rowUncoveredIdx, i)
    end

    @inbounds for i in eachindex(colCovered)
        colCovered[i] ? push!(colCoveredIdx, i) : push!(colUncoveredIdx, i)
    end

    h = typemax(T)
    @inbounds for j in colUncoveredIdx, i in rowUncoveredIdx
        cost = costMat[i,j] + Δcol[j] + Δrow[i]
        if cost <= h
            if cost != h
                h = cost
                empty!(minLocations)
            end
            push!(minLocations, (i,j))
        end
    end

    for i in rowCoveredIdx
        Δrow[i] += h
    end

    for i in colUncoveredIdx
        Δcol[i] -= h
    end

    # mark new zeros in Zs and remove those elements that will no longer be zero
    # produce new zeros
    for loc in minLocations
        Zs[loc...] = Z
    end
    # reduce old zeros
    rows = rowvals(Zs)
    for c in colCoveredIdx, i in nzrange(Zs, c)
        r = rows[i]
        if rowCovered[r]
            if Zs[r,c] == STAR
                rowSTAR[r] = false
                delete!(row2colSTAR, r)
            end
            Zs[r,c] = NON
        end
    end

    dropzeros!(Zs)

    # "return to step 1."
    return stepNum = 1
end


#### NOTICE!!!
#### most functions in this file are taken with almost no modifications from https://github.com/Gnimuc/Hungarian.jl
#### credit for those functions go to the original contributors: Yupei Qi <qiyupei@gmail.com>, et al
#### --------------------------------------------------------------------------------------------------------

###### ORIGINAL LICENSE
# MIT License
#
# Copyright (c) 2017 Yupei Qi and other contributors
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


end # module
