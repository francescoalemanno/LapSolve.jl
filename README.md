# LapSolve
Linear Assignment Problem solver via the Kuhn-Munkres algorithm, original Kernel functions come from the excellent [Hungarian.jl](https://github.com/Gnimuc/Hungarian.jl), this package adds a solver for stiff problems, which are LAP problems with too many blocking costs where an optimal fully assigned solution is impossible, but an optimal partial assignment may be possible. This solver has been modified to be completely non-copying and non-mutating, for saving as much memory as possible (this modification has also been sent upstream [#15](https://github.com/Gnimuc/Hungarian.jl/pull/15), although in some cases Hungarian.jl will still copy).

[![Build Status](https://travis-ci.com/francescoalemanno/LapSolve.jl.svg?branch=master)](https://travis-ci.com/francescoalemanno/LapSolve.jl)

[![codecov](https://codecov.io/gh/francescoalemanno/LapSolve.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/francescoalemanno/LapSolve.jl)

[![Coverage Status](https://coveralls.io/repos/github/francescoalemanno/LapSolve.jl/badge.svg?branch=master)](https://coveralls.io/github/francescoalemanno/LapSolve.jl?branch=master)

### Example
```julia
julia> using LapSolve

julia> M=rand(1:100,3,3)
3×3 Array{Int64,2}:
 66  83  92
 17  92  65
 30  84  18

julia> solve_lap(M)
([2, 1, 3], 118)
```

### Example for Stiff Problem
```julia
julia> using LapSolve

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
