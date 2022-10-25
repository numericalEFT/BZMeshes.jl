module BaseMesh

using ..StaticArrays
using ..LinearAlgebra

using ..BaryCheb
using ..AbstractMeshes
using ..Model

export UniformMesh, BaryChebMesh, CenteredMesh, EdgedMesh, UMesh

struct UMesh{T,DIM} <: AbstractMesh{T,DIM}
    lattice::Matrix{T}
    inv_lattice::Matrix{T}
    cell_volume::T
    origin::SVector{DIM,T}
    size::NTuple{DIM,Int}
    shift::SVector{DIM,Rational}
end

UMesh(;
    br::Brillouin{T,DIM},
    origin::Real,
    size,
    shift::Real) where {T,DIM} = UMesh{T,DIM}(
    br.recip_lattice,
    br.inv_recip_lattice,
    br.recip_cell_volume,
    SVector{DIM,T}(br.recip_lattice * ones(T, DIM) .* origin),
    size,
    SVector{DIM,Rational}(shift .* ones(Int, DIM))
)

Base.length(mesh::UMesh) = prod(mesh.size)
Base.size(mesh::UMesh) = mesh.size
Base.size(mesh::UMesh, I) = mesh.size[I]

function Base.show(io::IO, mesh::UMesh)
    println("UMesh with $(length(mesh)) mesh points")
end

function AbstractMeshes.fractional_coordinates(mesh::UMesh{T,DIM}, I::Int) where {T,DIM}
    n = SVector{DIM,Int}(ind2inds(mesh.size, I))
    return (n .- 1 .+ mesh.shift) ./ mesh.size
end

function AbstractMeshes.fractional_coordinates(mesh::UMesh{T,DIM}, x::AbstractVector) where {T,DIM}
    displacement = SVector{DIM,T}(x)
    return (mesh.inv_lattice * displacement)
end

function Base.getindex(mesh::UMesh{T,DIM}, inds...) where {T,DIM}
    n = SVector{DIM,Int}(inds)
    return mesh.origin + mesh.lattice * ((n .- 1 .+ mesh.shift) ./ mesh.size)
end

function Base.getindex(mesh::UMesh, I::Int)
    return Base.getindex(mesh, _ind2inds(mesh.size, I)...)
end

function AbstractMeshes.locate(mesh::UMesh{T,DIM}, x) where {T,DIM}
    # find index of nearest grid point to the point
    svx = SVector{DIM,T}(x)
    inds = fractional_coordinates(mesh, svx - mesh.origin) .* mesh.size .+ 1.5 .- mesh.shift .+ 2 .* eps.(T.(mesh.size))
    indexall = 1
    # println((mesh.invlatvec * displacement))
    # println(inds)
    factor = 1
    indexall += (_indfloor(inds[1], mesh.size[1]; edgeshift=0) - 1) * factor
    for i in 2:DIM
        factor *= mesh.size[i-1]
        indexall += (_indfloor(inds[i], mesh.size[i]; edgeshift=0) - 1) * factor
    end

    return indexall
end
AbstractMeshes.volume(mesh::UMesh) = mesh.volume
AbstractMeshes.volume(mesh::UMesh, i) = mesh.volume / length(mesh)


# in spglib, grid_address runs from 1-ceil(N/2) to N-ceil(N/2)
# thus -1:2 for N=4 and -2:2 for N=5

function spglib_grid_address_to_index(mesh::UMesh{T,DIM}, ga) where {T,DIM}
    inds = ga[1:DIM] # if length(x)==3 but DIM==2, take first two
    fcoords = (inds .+ mesh.shift) ./ mesh.size
    # shift fcoords 
    fcoords = [(fcoords[i] < 0) ? (fcoords[i] + 1) : fcoords[i] for i in 1:DIM]
    x = mesh.origin + mesh.lattice * fcoords
    return locate(mesh, x)
end
#####################################
# LEGACY CODE BELOW
#####################################

abstract type EqualLengthMesh{DIM,N} <: AbstractMesh{Float64,DIM} end

abstract type MeshType end
struct CenteredMesh <: MeshType end # Monkhorst-Pack mesh, take center points instead of left-bottom
struct EdgedMesh <: MeshType end # Γ=(0,0) centered, take (0,0) as mesh point

# TODO: support (N1, N2, N3)
struct UniformMesh{DIM,N,MT,DIMSQ} <: EqualLengthMesh{DIM,N}
    # DIMSQ should be provided explicitly in SMatrix parameters
    # otherwise it cause type instability
    origin::SVector{DIM,Float64}
    latvec::SMatrix{DIM,DIM,Float64,DIMSQ}
    invlatvec::SMatrix{DIM,DIM,Float64,DIMSQ}
    # dims could be here as field element

    function UniformMesh{DIM,N,MT}(origin, latvec) where {DIM,N,MT<:MeshType}
        return new{DIM,N,MT,DIM^2}(origin, latvec, inv(latvec))
    end
end

function UniformMesh{DIM,N}(origin, latvec) where {DIM,N}
    # Centered Mesh is the default
    return UniformMesh{DIM,N,CenteredMesh}(origin, latvec)
end

Base.length(mesh::UniformMesh{DIM,N}) where {DIM,N} = N^DIM
Base.size(mesh::UniformMesh{DIM,N}) where {DIM,N} = NTuple{DIM,Int}(ones(Int, DIM) .* N)

function Base.show(io::IO, mesh::UniformMesh)
    println("UniformMesh Grid:")
    for (pi, p) in enumerate(mesh)
        println(p)
    end
end

# Base.view(mesh::UniformMesh, i::Int) = Base.view(mesh.mesh, :, i)

# # set is not allowed for meshs
meshshift(::Type) = error("not implement!")
meshshift(::Type{<:CenteredMesh}) = 0.5
meshshift(::Type{<:EdgedMesh}) = 0.0

function Base.getindex(mesh::UniformMesh{DIM,N,MT}, inds...) where {DIM,N,MT}
    # pos = Vector(mesh.origin)
    # for (ni, n) in enumerate(inds)
    #     pos = pos .+ mesh.latvec[:, ni] .* (n - 1 + meshshift(MT)) ./ (N)
    # end
    # return pos
    n = SVector{DIM,Int}(inds)
    return mesh.origin + mesh.latvec * ((n .- 1 .+ meshshift(MT)) ./ N)
end

function Base.getindex(mesh::UniformMesh{DIM,N,MT}, i::Int) where {DIM,N,MT}
    return Base.getindex(mesh, _ind2inds(mesh, i)...)
end
Base.firstindex(mesh::UniformMesh) = 1
Base.lastindex(mesh::UniformMesh) = length(mesh)
# # iterator
Base.iterate(mesh::UniformMesh) = (mesh[1], 1)
Base.iterate(mesh::UniformMesh, state) = (state >= length(mesh)) ? nothing : (mesh[state+1], state + 1)

# _ind2inds(i::Int, N::Int, DIM::Int) = digits(i - 1, base=N, pad=DIM) .+ 1

# function _inds2ind(inds, N::Int)
#     indexall = 1
#     for i in 1:length(inds)
#         indexall += (inds[i] - 1) * N^(i - 1)
#     end
#     return indexall
# end

@generated function _inds2ind(umesh::EqualLengthMesh{DIM,N}, I) where {DIM,N}
    ex = :(I[DIM] - 1)
    for i = (DIM-1):-1:1
        ex = :(I[$i] - 1 + N * $ex)
    end
    return :($ex + 1)
end

@generated function _ind2inds(umesh::EqualLengthMesh{DIM,N}, I::Int) where {DIM,N}
    inds, quotient = :((I - 1) % N + 1), :((I - 1) ÷ N)
    for i = 2:DIM-1
        inds, quotient = :($inds..., $quotient % N + 1), :($quotient ÷ N)
    end
    inds = :($inds..., $quotient + 1)
    return :(SVector{DIM,Int}($inds))
end

function _indfloor(x, N; edgeshift=1)

    # edgeshift = 1 by default in floor function so that end point return N-1
    # edgeshift = 0 in locate function
    if x < 1
        return 1
    elseif x >= N
        return N - edgeshift
    else
        return floor(Int, x)
    end
end

function Base.floor(mesh::UniformMesh{DIM,N}, x) where {DIM,N}
    # find index of nearest grid point to the point
    displacement = SVector{DIM,Float64}(x) - mesh.origin
    # println(displacement)
    inds = (mesh.invlatvec * displacement) .* (N) .+ 0.5 .+ 2 * eps(1.0 * N)
    indexall = 1
    # println((mesh.invlatvec * displacement))
    # println(inds)
    for i in 1:DIM
        indexall += (_indfloor(inds[i], N) - 1) * N^(i - 1)
    end

    return indexall
end

function AbstractMeshes.locate(mesh::UniformMesh{DIM,N,MT}, x) where {DIM,N,MT}
    # find index of nearest grid point to the point
    displacement = SVector{DIM,Float64}(x) - mesh.origin
    # println(displacement)
    inds = (mesh.invlatvec * displacement) .* (N) .+ 1.5 .- meshshift(MT) .+ 2 * eps(1.0 * N)
    indexall = 1
    # println((mesh.invlatvec * displacement))
    # println(inds)
    for i in 1:DIM
        indexall += (_indfloor(inds[i], N; edgeshift=0) - 1) * N^(i - 1)
    end
    return indexall
end

AbstractMeshes.volume(mesh::UniformMesh) = abs(det(mesh.latvec))
function AbstractMeshes.volume(mesh::UniformMesh{DIM,N,MT}, i) where {DIM,N,MT}
    volume(MT, mesh, i)
end

function AbstractMeshes.volume(::Type{<:CenteredMesh}, mesh, i)
    # for uniform centered mesh, all mesh points share the same volume
    return abs(det(mesh.latvec)) / length(mesh)
end

function AbstractMeshes.volume(::Type{<:EdgedMesh}, mesh::UniformMesh{DIM,N,MT}, i) where {DIM,N,MT}
    inds = _ind2inds(mesh, i)
    n1, nend = count(i -> i == 1, inds), count(i -> i == N, inds)
    cellarea = 2^(DIM - n1 - nend) * 3^nend / 2^DIM
    return cellarea / length(mesh) * volume(mesh)
end

function interp(data, mesh::UniformMesh{DIM,N}, x) where {DIM,N}
    error("Not implemented!")
end

function interp(data, mesh::UniformMesh{2,N,MT}, x) where {N,MT}
    # find floor index and normalized x y
    displacement = SVector{2,Float64}(x) - mesh.origin
    xy = (mesh.invlatvec * displacement) .* N .+ (1 - meshshift(MT)) .* (1 + 2 * eps(N * 1.0))
    xi, yi = _indfloor(xy[1], N), _indfloor(xy[2], N)

    return linear2D(data, xi, yi, xy..., N)
end

@inline function linear2D(data, xi, yi, x, y, N)
    # accept data, floored index, normalized x and y, return linear interp
    # (xi, yi) should be [(1, 1) - size(data)], x and y normalized to the same scale as xi and yi
    xd, yd = x - xi, y - yi

    c0 = data[xi+(yi-1)*N] * (1 - xd) + data[xi+1+(yi-1)*N] * xd
    c1 = data[xi+(yi)*N] * (1 - xd) + data[xi+1+(yi)*N] * xd

    return c0 * (1 - yd) + c1 * yd
end

function interp(data, mesh::UniformMesh{3,N,MT}, x) where {T,N,MT}
    # find floor index and normalized x y z
    displacement = SVector{3,Float64}(x) - mesh.origin
    xyz = (mesh.invlatvec * displacement) .* N .+ (1 - meshshift(MT)) .* (1 + 4 * eps(N * 1.0))
    xi, yi, zi = _indfloor(xyz[1], N), _indfloor(xyz[2], N), _indfloor(xyz[3], N)

    return linear3D(data, xi, yi, zi, xyz..., N)
end

@inline function linear3D(data, xi, yi, zi, x, y, z, N) where {T}
    xd, yd, zd = x - xi, y - yi, z - zi

    c00 = data[xi+(yi-1)*N+(zi-1)*N^2] * (1 - xd) + data[xi+(yi-1)*N+(zi-1)*N^2] * xd
    c01 = data[xi+(yi-1)*N+(zi)*N^2] * (1 - xd) + data[xi+(yi-1)*N+(zi)*N^2] * xd
    c10 = data[xi+(yi)*N+(zi-1)*N^2] * (1 - xd) + data[xi+(yi)*N+(zi-1)*N^2] * xd
    c11 = data[xi+(yi)*N+(zi)*N^2] * (1 - xd) + data[xi+(yi)*N+(zi)*N^2] * xd

    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd

    return c0 * (1 - zd) + c1 * zd
end

function integrate(data, mesh::UniformMesh{DIM,N,CenteredMesh}) where {DIM,N}
    area = abs(det(mesh.latvec))
    avg = sum(data) / length(data)
    return avg * area
end

function integrate(data, mesh::UniformMesh{DIM,N,EdgedMesh}) where {DIM,N}
    area = abs(det(mesh.latvec))
    avg = 0.0
    for i in 1:length(data)
        inds = _ind2inds(mesh, i)
        n1, nend = count(i -> i == 1, inds), count(i -> i == N, inds)
        avg += 2^(DIM - n1 - nend) * 3^nend / 2^DIM * data[i]
    end
    avg = avg / length(data)
    return avg * area
end

struct BaryChebMesh{DIM,N,DIMSQ} <: EqualLengthMesh{DIM,N}
    origin::SVector{DIM,Float64}
    latvec::SMatrix{DIM,DIM,Float64,DIMSQ}
    invlatvec::SMatrix{DIM,DIM,Float64,DIMSQ}
    barycheb::BaryCheb1D{N} # grid in barycheb is [-1, 1]

    function BaryChebMesh{DIM,N}(origin, latvec, barycheb::BaryCheb1D{N}) where {DIM,N}
        return new{DIM,N,DIM^2}(origin, latvec, inv(latvec), barycheb)
    end
end

function BaryChebMesh(origin, latvec, DIM, N)
    barycheb = BaryCheb1D(N)
    return BaryChebMesh{DIM,N}(origin, latvec, barycheb)
end

function BaryChebMesh(origin, latvec, bcmesh::BaryChebMesh{DIM,N}) where {DIM,N}
    # this constructor borrows barycheb from generated bcmesh
    barycheb = bcmesh.barycheb
    return BaryChebMesh{DIM,N}(origin, latvec, barycheb)
end

Base.length(mesh::BaryChebMesh{DIM,N}) where {DIM,N} = N^DIM
Base.size(mesh::BaryChebMesh{DIM,N}) where {DIM,N} = NTuple{DIM,Int}(ones(Int, DIM) .* N)
function Base.show(io::IO, mesh::BaryChebMesh)
    println("BaryChebMesh Grid:")
    for (pi, p) in enumerate(mesh)
        println(p)
    end
end
function Base.getindex(mesh::BaryChebMesh{DIM,N}, inds...) where {DIM,N}
    # pos = Vector(mesh.origin)
    # for (ni, n) in enumerate(inds)
    #     pos = pos .+ mesh.latvec[:, ni] .* (mesh.barycheb[n] + 1.0) ./ 2.0
    # end
    # return pos
    a = SVector{DIM,Float64}(mesh.barycheb[n] for n in inds)
    return mesh.origin + mesh.latvec * (a .+ 1.0) ./ 2.0
end
function Base.getindex(mesh::BaryChebMesh{DIM,N}, i::Int) where {DIM,N}
    return Base.getindex(mesh, _ind2inds(mesh, i)...)
end

Base.firstindex(mesh::BaryChebMesh) = 1
Base.lastindex(mesh::BaryChebMesh) = length(mesh)
# # iterator
Base.iterate(mesh::BaryChebMesh) = (mesh[1], 1)
Base.iterate(mesh::BaryChebMesh, state) = (state >= length(mesh)) ? nothing : (mesh[state+1], state + 1)

function interp(data, mesh::BaryChebMesh{DIM,N}, x) where {DIM,N}
    # translate into dimensionless
    displacement = SVector{DIM,Float64}(x) - mesh.origin
    xs = (mesh.invlatvec * displacement) .* 2.0 .- 1.0
    return interpND(data, mesh.barycheb, xs)
end

function integrate(data, mesh::BaryChebMesh{DIM,N}) where {DIM,N}
    area = abs(det(mesh.latvec))
    return integrateND(data, mesh.barycheb, DIM) * area / 2^DIM
end

function locate1d(g::BaryCheb1D, x)
    grid = g.x
    if x <= grid[1]
        return 1
    elseif x >= grid[end]
        if length(grid) != 1
            return length(grid)
        end
    end

    i2 = searchsortedfirst(grid, x)
    i1 = i2 - 1
    return abs(grid[i1] - x) < abs(grid[i2] - x) ? i1 : i2
end

function volume1d(g::BaryCheb1D, i)
    grid = g.x
    if i != 1 && i != length(grid)
        return (grid[i+1] - grid[i-1]) / 2
    elseif i == 1
        return (grid[i+1] + grid[i]) / 2 - (-1)
    else
        return 1 - (grid[i] + grid[i-1]) / 2
    end
end

function AbstractMeshes.locate(mesh::BaryChebMesh{DIM,N}, x) where {DIM,N}
    displacement = SVector{DIM,Float64}(x) - mesh.origin
    xs = (mesh.invlatvec * displacement) .* 2.0 .- 1.0
    inds = [locate1d(mesh.barycheb, xs[i]) for i in 1:DIM]
    return _inds2ind(mesh, inds)
end

AbstractMeshes.volume(mesh::BaryChebMesh) = abs(det(mesh.latvec))
function AbstractMeshes.volume(mesh::BaryChebMesh{DIM,N}, i) where {DIM,N}
    inds = _ind2inds(mesh, i)
    return reduce(*, volume1d(mesh.barycheb, inds[j]) for j in 1:DIM) * volume(mesh) / 2^DIM
end

end

