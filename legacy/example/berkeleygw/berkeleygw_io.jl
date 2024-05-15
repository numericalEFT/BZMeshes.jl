using JLD2
using BrillouinZoneMeshes
using UnPack
using StaticArrays

struct GVectors{DIM} <: AbstractMeshes.AbstractMesh{Int,DIM}
    gmin::SVector{DIM,Int}
    # gmax::SVector{DIM,Int} gmax = gmin .+ size .- 1
    size::NTuple{DIM,Int}
end

function GVectors(gmin, gmax)
    DIM = length(gmin)
    @assert DIM == length(gmax)
    size = tuple((gmax .- gmin .+ 1)...)
    return GVectors{DIM}(gmin, size)
end

Base.getindex(gvs::GVectors, inds...) = gvs.gmin .+ inds .- 1
Base.getindex(gvs::GVectors, I::Int) = gvs[AbstractMeshes._ind2inds(size(gvs), I)...]
Base.show(io::IO, mesh::GVectors) = print(io, "GVectors of $(mesh.size)")

_locate(gvs::GVectors, gv) = AbstractMeshes._inds2ind(size(gvs), (gv .- gvs.gmin .+ 1))
function locate(gvs::GVectors, gv)
    # return 0 if not in gvs
    # return index otherwise
    I = _locate(gvs, gv)
    if 0 < I <= length(gvs)
        return I
    else
        return 0
    end
end

function load_cell(wfn)
    # load cell info from wfn.h5
    @unpack celvol, recvol, alat, blat, nat, avec, bvec, adot, bdot, atyp, apos = wfn["mf_header/crystal"]
    println(alat * blat)
    println(celvol * recvol)
    println(avec)
    # cell = Cells.Cell(;lattice=)
end

function load_rbzmesh(wfn)
    # load cell
    bgw_lattice = wfn[""]

    bgw_kgrid = wfn["mf_header/kpoints/kgrid"]
    bgw_shift = wfn["mf_header/kpoints/shift"]

end

if abspath(PROGRAM_FILE) == @__FILE__

    wfn = jldopen("wfn.h5", "r")
    cell = load_cell(wfn)
    # rbzmesh = load_rbzmesh(wfn)
    # println(rbzmesh)

end