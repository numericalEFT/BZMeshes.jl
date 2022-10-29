using BrillouinZoneMeshes.PointSymmetry: spglib_spacegroup_number, spglib_standardize_cell

@testset "Brillouin" begin
    Brillouin = BrillouinZoneMeshes.Model.Brillouin

    # square lattice
    DIM = 2
    lattice = Matrix([1.0 0; 0 1]')
    br = Brillouin(lattice=lattice)
    @test br.inv_lattice .* 2π ≈ br.recip_lattice'
    @test br.unit_cell_volume ≈ abs(det(lattice))
    @test br.recip_cell_volume ≈ 1 / abs(det(lattice)) * (2π)^DIM

    # triagular lattice
    DIM = 2
    lattice = Matrix([2.0 0; 1 sqrt(3)]')
    br = Brillouin(lattice=lattice)
    @test br.inv_lattice .* 2π ≈ br.recip_lattice'
    @test br.unit_cell_volume ≈ abs(det(lattice))
    @test br.recip_cell_volume ≈ 1 / abs(det(lattice)) * (2π)^DIM

    # 3d testing lattice
    DIM = 3
    lattice = Matrix([2.0 0 0; 1 sqrt(3) 0; 7 11 19]')
    br = Brillouin(lattice=lattice)
    @test br.inv_lattice .* 2π ≈ br.recip_lattice'
    @test br.unit_cell_volume ≈ abs(det(lattice))
    @test br.recip_cell_volume ≈ 1 / abs(det(lattice)) * (2π)^DIM
end

@testset "Standard Brillouin" begin
    a = 10.3
    Si = 1
    Ge = 2

    # silicon with Cartesian x coordinates flipped
    lattice = a / 2 * [[0 -1 -1.0]; [1 0 1.0]; [1 1 0.0]]
    atoms = [Si, Si]
    positions = [ones(3) / 8, -ones(3) / 8]
    model = BrillouinZoneMeshes.Model.standard_brillouin(lattice=lattice, atoms=atoms, positions=positions, primitive=false)
    println(model)
    display(model)

    @test spglib_spacegroup_number(model) == 227
    @test model.lattice ≈ a * I(3)

    # Zincblende structure with different lattice vectors
    lattice = a / 2 * [[0 1 1.0]; [-1 0 1.0]; [-1 1 0.0]]
    atoms = [Si, Ge]
    positions = [[-1, 1, 1] / 8, -[-1, 1, 1] / 8]
    model = BrillouinZoneMeshes.Model.standard_brillouin(lattice=lattice, atoms=atoms, positions=positions, primitive=false)
    @test spglib_spacegroup_number(model) == 216
    @test model.lattice ≈ a * I(3)

    # Two-dimensional example
    lattice = [[1.0 0.0]; [0.0 1.0]]
    atoms = [Si,]
    positions = [[0.0, 0.0],]
    model = BrillouinZoneMeshes.Model.standard_brillouin(lattice=lattice, atoms=atoms, positions=positions, primitive=true)
    @test model.lattice ≈ lattice
end
