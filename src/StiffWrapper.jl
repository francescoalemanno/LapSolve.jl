
struct StiffWrapper{T,S} <: AbstractMatrix{T}
    parent::S
    max::T
    min::T
    penalty::T
    translation::T
    @inline function StiffWrapper{T,S}(A::S,penalty) where {T,S}
        maxV=typemin(T)
        minV=typemax(T)
        for i in eachindex(A)
            v=A[i]
            isinf(v) && continue
            maxV=max(maxV,v)
            minV=min(minV,v)
        end
        if isinf(maxV) || isinf(minV)
            maxV=minV=zero(T)
        end
        trasl=-minV
        if iszero(maxV+trasl)
            trasl+=one(T) #regularize illposed stiff problems
        end
        new(A,maxV,minV,penalty,trasl)
    end
end
@inline StiffWrapper(A,penalty)=StiffWrapper{eltype(A),typeof(A)}(A,penalty)

import Base.size, Base.getindex, Base.IndexStyle

@inline IndexStyle(::StiffWrapper) = IndexCartesian()
@inline function size_parent(A::StiffWrapper)
    size(A.parent)
end
@inline function size(A::StiffWrapper)
    n,m=size_parent(A)
    (n+m,n+m)
end
@inline function getindex(A::StiffWrapper{T},i::Int,j::Int) where {T}
    n,m=size_parent(A)
    @fastmath begin
    if i<=n && j<=m
        return A.parent[i,j]+A.translation
    end
    if i==j+n || j==i+m
        return (A.max+A.translation)*A.penalty
    end
    if i>n && j>m
        v=A.parent[j-m,i-n]+A.translation
        if !isinf(v)
            return A.min+A.translation
        end
    end
    return typemax(T)
end
end
