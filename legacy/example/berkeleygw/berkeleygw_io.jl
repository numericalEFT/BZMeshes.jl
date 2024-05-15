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
    lattice = alat .* avec
    # currently only for simple case
    # atoms = Vector{Int64}(atyp)
    atoms = Vector{Int64}([i for i in 1:length(atyp)])
    positions = [apos[:, i] for i in 1:nat]
    cell = Cells.Cell(; lattice=Matrix(lattice), atoms=atoms, positions=positions)
    return cell
end

function compute_meshmap(bzmesh, rk, mtrx)
    nrk = size(rk)[2]
    irreducible_indices = zeros(Int, nrk)
    for i in 1:nrk
        k = rk[:, i]
        irreducible_indices[i] = AbstractMeshes.locate(bzmesh, AbstractMeshes.frac_to_cart(bzmesh, k))
        @assert bzmesh[AbstractMeshes.FracCoords, irreducible_indices[i]] ≈ k
    end
    # constructing map
    map = zeros(Int, length(bzmesh))
    for i in 1:length(bzmesh)
        k_frac = Vector(bzmesh[AbstractMeshes.FracCoords, i])
        for m in 1:3
            if k_frac[m] >= 0.5
                k_frac[m] -= 1
            end
        end
        k = AbstractMeshes.frac_to_cart(bzmesh, k_frac)
        for j in 1:size(mtrx)[3]
            symat = mtrx[:, :, j]
            newk = symat * k
            newk_frac = AbstractMeshes.cart_to_frac(bzmesh, newk)
            for m in 1:3
                if newk_frac[m] >= 0.5
                    newk_frac[m] -= 1
                elseif newk_frac[m] < -0.5
                    newk_frac[m] += 1
                end
            end
            newk_frac .+= 0.5
            newk = AbstractMeshes.frac_to_cart(bzmesh, newk_frac)
            for m in 1:nrk
                frac_rk = rk[:, m]
                cart_rk = AbstractMeshes.frac_to_cart(bzmesh, frac_rk)
                if isapprox(cart_rk, newk, atol=1e-6)
                    # println("k=$k, mat=$(symat), newk=$newk, rk=$(cart_rk)")
                    println("k=$k_frac, mat=$(symat), newk=$newk_frac, rk=$(frac_rk)")
                    # @assert (map[i] == 0 || map[i] == m) "$k reduced to two points:$(map[i]) and $m !"
                    map[i] = m
                end
            end
        end
        # @assert map[i] != 0 "found no rk for i=$i, k=$k !"
    end
    return
end

function load_rbzmesh(wfn)
    # load cell
    cell = load_cell(wfn)
    @unpack ntran, cell_symmetry, mtrx, tnp = wfn["mf_header/symmetry"]
    @unpack nspin, nspinor, nrk, mnband, ngkmax, ecutwfc, kgrid, shift, ngk, ifmin, ifmax, w, rk, el, occ, = wfn["mf_header/kpoints"]

    isshift = [(s == 0.5) ? true : false for s in shift] # only apply for 0 and 1/2 shift
    # BGW set origin at 0
    bzmesh = BZMeshes.UniformBZMesh(cell=cell, origin=0.0, size=kgrid, shift=isshift)

    # construct symmetry map
    meshmap = compute_meshmap(bzmesh, rk, mtrx)
end

if abspath(PROGRAM_FILE) == @__FILE__

    wfn = jldopen("wfn.h5", "r")
    cell = load_cell(wfn)
    # rbzmesh = load_rbzmesh(wfn)
    # println(rbzmesh)

end