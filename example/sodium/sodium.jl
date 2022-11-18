using DFTK
using Unitful
using UnitfulAtomic
using Plots
# using PlotlyJS
using Printf
using PyCall
using BrillouinZoneMeshes

include("interpolation.jl")

# Setup a sodium lattice, see the following link for more information
# http://lampx.tugraz.at/~hadley/ss1/crystalstructure/structures/bcc/bcc.php
a = 0.429u"nm"

# each column is a lattice vector
lattice = a / 2 * [[1.0 -1.0 1.0]; [1.0 1.0 -1.0]; [-1.0 1.0 1.0]]
pseudo_file = joinpath(@__DIR__, "Na_PBE.upf")
Na = ElementPsp(:Na, psp=load_psp(pseudo_file))
atoms = [Na,]
positions = [[0.0, 0.0, 0.0],]

smearing = DFTK.Smearing.Gaussian()
temperature = 0.02u"Ry"

model = model_PBE(lattice, atoms, positions; temperature, smearing)

Ecut = 80.0u"Ry"
kgrid = [16, 16, 16]

basis = PlaneWaveBasis(model; Ecut, kgrid=kgrid)

scfres = self_consistent_field(basis, tol=1e-8);

# plot
plot_dos(scfres, unit=u"eV", xlims=[-5, 20], ylims=[0, 40]) # in hartree unit
plot_bandstructure(scfres, kline_density=10, unit=u"eV", ρ=scfres.ρ)

# calculate band width
_, data = compute_bands(basis, [[0.0, 0.0, 0.0],], ρ=scfres.ρ);

bandwidth = uconvert.(u"eV", (data[1][5] - scfres.εF) * u"hartree")
println(bandwidth)

# # interpolation of bands
# interp = pyimport("scipy.interpolate")
# np = pyimport("numpy")
# griddata = interp.griddata

# build BZ mesh and maps
cell = Cells.Cell(; lattice=austrip.(lattice), atoms=[1,], positions)
bzmesh = BZMeshes.Monkhorst_Pack(cell=cell, size=kgrid, shift=[false, false, false])
meshmap = MeshMaps.MeshMap(bzmesh)

println(length(meshmap.irreducible_indices))
for idx in meshmap.irreducible_indices
    k = inv_lattice_vector(bzmesh) * bzmesh[idx]
    println(k)
end

# band structure with BZ mesh and maps
band = Dict{Int,Float64}()
bandix = 5 # the fifth band is the valence band
for (ki, kpoints) in enumerate(basis.kpoints)
    idx = AbstractMeshes.locate(bzmesh, lattice_vector(bzmesh) * basis.kpoints[ki].coordinate)
    band[meshmap.map[idx]] = scfres.eigenvalues[ki][bandix] - scfres.εF
end

# band fourier transform
bandarray = zeros(kgrid...);

for kx in 1:kgrid[1]
    for ky in 1:kgrid[2]
        for kz in 1:kgrid[3]
            idx = AbstractMeshes._inds2ind(tuple(kgrid...), [kx, ky, kz])
            # idx = AbstractMeshes.locate(bzmesh, bzmesh[kx, ky, kz])
            bandarray[kx, ky, kz] = band[meshmap.map[idx]]
        end
    end
end

Kfine = 100
band_interpolate = fourier_interpolate(bandarray, Kfine);
Plots.plot()
Plots.plot!((0:kgrid[1]-1) / kgrid[1], bandarray[:, 1, 1], label="band", markershape=:circle)
Plots.plot!((0:Kfine-1) / Kfine, band_interpolate[:, 1, 1], label="band interpolation")

band_interpolate = fourier_interpolate(bandarray[:, 8, 8], Kfine);
Plots.plot()
Plots.plot!((0:kgrid[1]-1) / kgrid[1], bandarray[:, 8, 8], label="band", markershape=:circle)
Plots.plot!((0:Kfine-1) / Kfine, band_interpolate, label="band interpolation")