
struct StiffWrapper{T,S} <: AbstractMatrix{T}
    parent::S
    max::T
    min::T
    penalty
    function StiffWrapper{T,S}(A::S,penalty) where {T,S}
        maxV=typemin(T)
        minV=typemax(T)
        for i in eachindex(A)
            v=A[i]
            isinf(v)&&continue
            maxV=max(maxV,v)
            minV=min(minV,v)
        end
        if isinf(maxV) || isinf(minV)
            error("Problem is unsolvable")
        end
        new(A,maxV,minV,convert(T,penalty))
    end
end
StiffWrapper(A,penalty)=StiffWrapper{eltype(A),typeof(A)}(A,penalty)

import Base.size, Base.getindex, Base.IndexStyle

IndexStyle(::Type{StiffWrapper}) = IndexCartesian()
function size_parent(A::StiffWrapper)
    size(A.parent)
end
function size(A::StiffWrapper)
    n,m=size_parent(A)
    (n+m,n+m)
end
function getindex(A::StiffWrapper{T},i::Int,j::Int) where {T}
    n,m=size_parent(A)
    if i<=n && j<=m
        return A.parent[i,j]
    end
    if i==j+n || j==i+m
        return A.max*A.penalty
    end
    if i>n && j>m
        v=A.parent[j-m,i-n]
        if !isinf(v)
            return A.min
        end
    end
    return typemax(T)
end
