
include("berkeleygw_io.jl")

using Test

@testset "Cell" begin

    wfn = jldopen("wfn.h5", "r")

    @testset "Test BerkeleyGW conventions" begin
        @unpack celvol, recvol, alat, blat, nat, avec, bvec, adot, bdot, atyp, apos = wfn["mf_header/crystal"]
        @test (alat * blat) ≈ 2π
        @test (celvol * recvol) ≈ (2π)^3
        println(alat)
    end

    # wfn = jldopen("wfn.h5", "r")
    # cell = load_cell(wfn)

end
