module DFTGreen

using DFTK
using Unitful
using UnitfulAtomic
using Interpolations
using BrillouinZoneMeshes
using JLD2

using BrillouinZoneMeshes.StaticArrays
using BrillouinZoneMeshes.LinearAlgebra

include("interpolation.jl")
using .FFTInterp

export fourier_interpolate
export calc_scf, save_scfres, load_scfres
export GVectors, locate
export GreenInterpolator, green

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
# compute scf, default sodium
function calc_scf(;
    lattice=(0.429u"nm") / 2 * [[1.0 -1.0 1.0]; [1.0 1.0 -1.0]; [-1.0 1.0 1.0]],
    pseudo_file=joinpath(@__DIR__, "Na_PBE.upf"),
    atoms=[ElementPsp(:Na, psp=load_psp(pseudo_file)),],
    positions=[[0.0, 0.0, 0.0],],
    smearing=DFTK.Smearing.Gaussian(),
    temperature=0.02u"Ry",
    Ecut=80.0u"Ry",
    kgrid=[16, 16, 16],
    tol=1e-10)

    model = model_PBE(lattice, atoms, positions; temperature, smearing)
    basis = PlaneWaveBasis(model; Ecut, kgrid=kgrid)
    scfres = self_consistent_field(basis, tol=tol)
    return scfres
end

function collect_gvector(scfres)
    gvectors = Vector{eltype(scfres.basis.G_vectors)}([])
    for (ik, k) in enumerate(scfres.basis.kpoints)
        for (ig, gv) in enumerate(k.G_vectors)
            if gv ∉ gvectors
                push!(gvectors, gv)
            end
        end
    end
    return gvectors
end

function gvector_minmax(gvectors)
    DIM = 3
    gmin = zeros(Int, DIM)
    gmax = zeros(Int, DIM)
    for (ig, gv) in enumerate(gvectors)
        for j in 1:DIM
            if gv[j] < gmin[j]
                gmin[j] = gv[j]
            end
            if gv[j] > gmax[j]
                gmax[j] = gv[j]
            end
        end
    end
    return gmin, gmax
end

function greenfromscf(scfres, m, n, kidx, τ, β)
    result = 0.0

    for band in 1:scfres.n_bands_converge
        ε = scfres.eigenvalues[kidx][band] - scfres.εF
        result += scfres.ψ[kidx][m, band] * conj(scfres.ψ[kidx][n, band]) * exp(-ε * τ) / (1 + exp(-ε * β))
    end
    return result
end

struct GreenInterpolator{DI,MT}
    # store interpolators for green's function and provide interpolation
    # interpolators are for k-points, thus for each G-vector and band an interpolator is needed
    # store ψ and ε, green's function is computed
    beta::Float64
    n_bands::Int
    dispersions::Vector{DI}
    ψ::Array{ComplexF64,3}
    gvectors::GVectors{3}
    rbzmesh::MT
end

function GreenInterpolator(scfres; n_bands=10)
    # extract parameters
    beta = scfres.εF / scfres.basis.model.temperature
    # n_bands = scfres.n_bands_converge
    n_bands = n_bands
    kgrid = scfres.basis.kgrid
    lattice = scfres.basis.model.lattice
    atoms = ones(Int, length(scfres.basis.model.atoms))
    positions = scfres.basis.model.positions

    # generate bzmesh and map
    cell = Cells.Cell(; lattice=Matrix(austrip.(lattice)), atoms=atoms, positions=positions)
    bzmesh = BZMeshes.Monkhorst_Pack(cell=cell, size=kgrid.kgrid_size, shift=[false, false, false])
    meshmap = MeshMaps.MeshMap(bzmesh)
    rbzmesh = MeshMaps.ReducedBZMesh(bzmesh, meshmap)

    kx = LinRange(-0.5, 0.5, kgrid.kgrid_size[1] + 1)
    ky = LinRange(-0.5, 0.5, kgrid.kgrid_size[2] + 1)
    kz = LinRange(-0.5, 0.5, kgrid.kgrid_size[3] + 1)


    gvfromscf = collect_gvector(scfres)
    gmin, gmax = gvector_minmax(gvfromscf)
    gvectors = GVectors(gmin, gmax)

    # collect data

    ψ = zeros(ComplexF64, (length(gvectors), n_bands, length(rbzmesh)))
    ϵ = zeros(Float64, (n_bands, length(rbzmesh)))
    for bandidx in 1:n_bands
        for ki in 1:length(scfres.basis.kpoints)
            k = scfres.basis.kpoints[ki]
            idx = AbstractMeshes.locate(rbzmesh, lattice_vector(rbzmesh.mesh) * k.coordinate)
            ϵ[bandidx, idx] = scfres.eigenvalues[ki][bandidx] - scfres.εF
            for gi in 1:length(k.G_vectors)
                g = k.G_vectors[gi]
                gidx = locate(gvectors, g)
                # @assert gvectors[gidx] ≈ g
                ψ[gidx, bandidx, idx] = scfres.ψ[ki][gi, bandidx]
            end
        end
    end

    # generate interpolator

    dispersions = []

    for bandidx in 1:n_bands
        # dispersion
        bandarray = zeros((kgrid.kgrid_size .+ 1)...)
        for ikx in 1:kgrid.kgrid_size[1]+1
            for iky in 1:kgrid.kgrid_size[2]+1
                for ikz in 1:kgrid.kgrid_size[3]+1
                    # idx = AbstractMeshes._inds2ind(tuple(kgrid.kgrid_size...), [kx, ky, kz])
                    idx = AbstractMeshes.locate(rbzmesh, lattice_vector(rbzmesh.mesh) * [kx[ikx], ky[iky], kz[ikz]])
                    bandarray[ikx, iky, ikz] = ϵ[bandidx, idx]
                end
            end
        end
        # itp = interpolate(bandarray, BSpline(Cubic(Line(OnGrid()))))
        itp = interpolate(bandarray, BSpline(Linear()))
        sitp = scale(itp, kx, ky, kz)
        push!(dispersions, sitp)
    end

    return GreenInterpolator{typeof(dispersions[1]),typeof(rbzmesh)}(beta, n_bands, dispersions,
        ψ, gvectors, rbzmesh)
end

function greenτ(gi::GreenInterpolator{DI}, gi1::Int, gi2::Int, k, τ;
    beta=gi.beta) where {DI}
    result = ComplexF64(0.0)
    # for band in 1:gi.n_bands
    for band in 1:gi.n_bands
        data1 = view(gi.ψ, gi1, band, :)
        ψ1 = AbstractMeshes.interp(data1, gi.rbzmesh, k)
        if gi1 == gi2
            ψ2 = ψ1
        else
            data2 = view(gi.ψ, gi2, band, :)
            ψ2 = AbstractMeshes.interp(data2, gi.rbzmesh, k)
        end
        sitp = gi.dispersions[band]
        fk = inv_lattice_vector(gi.rbzmesh.mesh) * k
        ε = sitp(fk[1], fk[2], fk[3])
        # println("ψ1=$ψ1, ψ2=$ψ2, ε=$ε")
        result += ψ1 * conj(ψ2) * exp(-ε * τ) / (1 + exp(-ε * beta))
    end
    return result
end

function greenτ(
    gi::GreenInterpolator{DI},
    gv1::AbstractVector{Int},
    gv2::AbstractVector{Int},
    k, τ;
    beta=gi.beta) where {DI}
    gi1, gi2 = locate(gi.gvectors, gv1), locate(gi.gvectors, gv2)
    if gi1 == 0 || gi2 == 0
        return ComplexF64(0.0)
    else
        return greenτ(gi, gi1, gi2, k, τ)
    end
end

function greenω(gi::GreenInterpolator{DI}, gi1::Int, gi2::Int, k, ω; δ=1e-8) where {DI}
    result = ComplexF64(0.0)
    for band in 1:gi.n_bands
        # for band in 5:5
        data1 = view(gi.ψ, gi1, band, :)
        ψ1 = AbstractMeshes.interp(data1, gi.rbzmesh, k)
        if gi1 == gi2
            ψ2 = ψ1
        else
            data2 = view(gi.ψ, gi2, band, :)
            ψ2 = AbstractMeshes.interp(data2, gi.rbzmesh, k)
        end
        sitp = gi.dispersions[band]
        fk = inv_lattice_vector(gi.rbzmesh.mesh) * k
        ε = sitp(fk[1], fk[2], fk[3])
        # println("ψ1=$ψ1, ψ2=$ψ2, ε=$ε")
        result += ψ1 * conj(ψ2) / (ω - ε + im * δ)
        # result += 1 / (ω - ε + im * δ) / length(gi.gvectors)
    end
    return result
end

function greenω(
    gi::GreenInterpolator{DI},
    gv1::AbstractVector{Int},
    gv2::AbstractVector{Int},
    k, ω; δ=1e-8) where {DI}
    gi1, gi2 = locate(gi.gvectors, gv1), locate(gi.gvectors, gv2)
    if gi1 == 0 || gi2 == 0
        return ComplexF64(0.0)
    else
        return greenω(gi, gi1, gi2, k, ω; δ)
    end
end

function greenfreeω(gi, gv1, gv2, k, ω; δ=1e-2)
    gi1, gi2 = locate(gi.gvectors, gv1), locate(gi.gvectors, gv2)
    if gi1 != gi2
        return ComplexF64(0.0)
    else
        K = k .+ lattice_vector(gi.rbzmesh.mesh) * gv1
        m = 1.0
        EF = 0.017217600764962135
        ε = dot(K, K) / 2 / m - EF
        return 1 / (ω - ε + im * δ)
    end
end

function greenfreeτ(gi, gi1::Int, gi2::Int, k, τ; beta=gi.beta)
    if gi1 != gi2
        return ComplexF64(0.0)
    else
        K = k .+ lattice_vector(gi.rbzmesh.mesh) * gi.gvectors[gi1]
        m = 0.5
        EF = 0.017217600764962135
        # EF = 0.0575495 # rs=8
        # EF = 0.230198 # rs=4
        # EF = 0.920792 # rs=2
        # EF = 3.68316855 # rs=1
        ε = dot(K, K) / 2 / m - EF
        return exp(-ε * τ) / (1 + exp(-ε * beta))
    end
end
function greenfreeτ(gi, gv1, gv2, k, τ; beta=gi.beta)
    gi1, gi2 = locate(gi.gvectors, gv1), locate(gi.gvectors, gv2)
    return greenfreeτ(gi, gi1, gi2, k, τ; beta=beta)
end

end

if abspath(PROGRAM_FILE) == @__FILE__

    using .DFTGreen
    scfres = DFTGreen.calc_scf(kgrid=[16, 16, 16])
    DFTGreen.save_scfres("./run/sodium.jld2", scfres)

end