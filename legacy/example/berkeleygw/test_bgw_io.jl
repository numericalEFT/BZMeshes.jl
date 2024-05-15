
include("berkeleygw_io.jl")

using Test, LinearAlgebra

@testset "Wave function" begin
    wfn = jldopen("wfn.h5", "r")

    @testset "Cell" begin

        @unpack celvol, recvol, alat, blat, nat, avec, bvec, adot, bdot, atyp, apos = wfn["mf_header/crystal"]

        @testset "Test BerkeleyGW conventions" begin
            @test celvol ≈ det(alat .* avec)
            @test recvol ≈ det(blat .* bvec)
            @test (alat * blat) ≈ 2π
            @test (celvol * recvol) ≈ (2π)^3

            @test nat == length(atyp) == size(apos)[2]
        end

        @testset "Test loaded cell" begin
            cell = load_cell(wfn)
            println(cell)
            @test cell.cell_volume ≈ celvol
            @test cell.recip_cell_volume ≈ recvol
            @test cell.lattice ≈ (alat .* avec)
            @test cell.recip_lattice ≈ (blat .* bvec)
        end

    end

    @testset "ReducedBZMesh" begin
        @unpack nspin, nspinor, nrk, mnband, ngkmax, ecutwfc, kgrid, shift, ngk, ifmin, ifmax, w, rk, el, occ, = wfn["mf_header/kpoints"]

        rbzmesh = load_rbzmesh(wfn)

    end

end
