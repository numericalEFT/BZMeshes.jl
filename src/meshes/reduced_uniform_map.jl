"""
Default logic to determine the symmetry operations to be used in the model.
"""
function default_symmetries(model::Model.Brillouin{T,DIM}
    ; tol_symmetry=PointSymmetry.SYMMETRY_TOLERANCE) where {T,DIM}

    lattice = zeros(T, 3, 3)
    for i in 1:DIM
        lattice[i, 1:DIM] = model.lattice[i, 1:DIM]
    end
    for i in DIM+1:3
        lattice[i, i] = 1
    end
    positions = [zeros(T, 3) for i in 1:length(model.atoms)]
    for ai in eachindex(model.atoms)
        positions[ai][1:DIM] = model.positions[ai][1:DIM]
    end
    atoms = model.atoms

    # Standard case from here on:
    if length(positions) != length(atoms)
        error("Length of atoms and positions vectors need to agree.")
    end

    return PointSymmetry.symmetry_operations(lattice, atoms, positions; tol_symmetry)
end


function _kcoords2ind(kcoord, kgrid_size, kshift)
    # kidx = [Int(kvec * kgrid_size[d] - kshift[d]) + kgrid_size[d] - 1 for (d, kvec) in enumerate(kcoord)]
    # inexact convert is not allowed with Int(), use floor(Int,) instead
    kidx = [floor(Int, (kvec + 1 / 2) * kgrid_size[d] - kshift[d]) for (d, kvec) in enumerate(kcoord)]
    kidx = [(kidx[d] + kgrid_size[d]) % kgrid_size[d] + 1 for d in 1:length(kidx)]
    klinearidx = AbstractMeshes._inds2ind(tuple(kgrid_size...), tuple(kidx...))
    return klinearidx
end

# WARNINING: Do not support Gamma_centered: origin=0  !!!
# Monkhorst-Pack: origin=-1/2, consistent with VASP
# to be consistent with DFTK: 
#  - N is even, VASP is the same as DFTK: shift=0 will include Gamma point, shift=1/2 will not
#  - N is odd, VASP is different as DFTK: shift=0 will not include Gamma point, shift=1/2 will
function _reduced_uniform_meshmap(model::Model.Brillouin{T,DIM}, symmetry::Bool=true;
    kgrid_size::Vector{Int}, kshift::Bool=false,
    tol_symmetry=PointSymmetry.SYMMETRY_TOLERANCE
) where {T,DIM}
    # Determine symmetry operations to use
    if symmetry
        symmetries = default_symmetries(model, tol_symmetry=tol_symmetry)
    else
        symmetries = [one(PointSymmetry.SymOp)]
    end
    @assert !isempty(symmetries)  # Identity has to be always present.
    _kgrid_size = ones(Int, 3)
    _kgrid_size[1:DIM] = kgrid_size[1:DIM]
    _kshift = kshift ? [1 // 2, 1 // 2, 1 // 2] : [0, 0, 0]

    kcoords, kweights, symmetries = PointSymmetry.bzmesh_ir_wedge(_kgrid_size, symmetries; kshift=_kshift)
    all_kcoords = PointSymmetry.unfold_kcoords(kcoords, symmetries)
    Nk = reduce(*, kgrid_size)
    @assert length(all_kcoords) == Nk
    kindices = []
    kmap = zeros(Int, Nk)
    inv_kmap = Dict{Int,Vector{Int}}()
    count = 0
    for kpoint in kcoords
        all_kpoint = PointSymmetry.unfold_kcoords([kpoint,], symmetries)
        k0ind = _kcoords2ind(kpoint, kgrid_size, _kshift)
        push!(kindices, k0ind)
        inv_kmap[k0ind] = []
        count += length(all_kpoint)
        for k in all_kpoint
            kind = _kcoords2ind(k, kgrid_size, _kshift)
            kmap[kind] = k0ind
            push!(inv_kmap[k0ind], kind)
        end
    end

    @assert count == Nk

    return MeshMaps.MeshMap(kindices, kmap, inv_kmap)
end

# in spglib, grid_address runs from 1-ceil(N/2) to N-ceil(N/2)
# thus -1:2 for N=4 and -2:2 for N=5

function spglib_grid_address_to_index(mesh::UMesh{T,DIM}, ga) where {T,DIM}
    inds = ga[1:DIM] # if length(x)==3 but DIM==2, take first two
    fcoords = (inds .+ mesh.shift) ./ mesh.size #fractional coordinates as defined in spglib
    # shift fcoords, nomalize to [0, 1)
    fcoords = [(fcoords[i] < 0) ? (fcoords[i] + 1) : fcoords[i] for i in 1:DIM]
    x = mesh.lattice * fcoords
    return locate(mesh, x)
end

function uniform_meshmap(mesh::BZMeshes.UniformBZMesh{T,DIM},
    symmetry::Bool=true;
    is_time_reversal::Bool=true,
    tol_symmetry=PointSymmetry.SYMMETRY_TOLERANCE
) where {T,DIM}
    # Determine symmetry operations to use
    if symmetry
        symmetries = default_symmetries(mesh.br, tol_symmetry=tol_symmetry)
    else
        symmetries = [one(PointSymmetry.SymOp)]
    end
    @assert mesh.mesh.inv_lattice * mesh.mesh.origin ≈ -ones(DIM) / 2 "Uniform BZ mesh MeshMap only supports origin=[-1/2, -1/2,...]"

    kgrid_size = mesh.mesh.size
    kshift = mesh.mesh.shift
    @assert !isempty(symmetries)  # Identity has to be always present.

    _kgrid_size = ones(Int, 3)
    _kgrid_size[1:DIM] .= kgrid_size[1:DIM]
    is_shift = [0, 0, 0]
    is_shift[1:DIM] .= Int.(kshift * 2)
    # kcoords_mp = [mesh.mesh.inv_lattice * mesh[k] for k in 1:length(mesh)]
    # println("sym before: ", length(symmetries))
    # return kcoords_mp

    #TODO: It is probably makes sense to use the symmetry operations to reduce the number of kpoints
    # Filter those symmetry operations that preserve the MP grid
    # kcoords_mp = kgrid_monkhorst_pack(kgrid_size; kshift)
    # symmetries = symmetries_preserving_kgrid(symmetries, kcoords_mp)

    # Ws = [symop.W for symop in symmetries]
    # _, mapping, grid = PointSymmetry.spglib_get_stabilized_reciprocal_mesh(
    #     kgrid_size, Ws; is_shift, is_time_reversal=false
    # )

    # lat, atoms, pos, mag_moments = PointSymmetry.spglib_standardize_cell(mesh.br, primitive=true)
    lat, atoms, pos, mag_moments = mesh.br.lattice, mesh.br.atoms, mesh.br.positions, []

    lat, pos = PointSymmetry._make3D(lat, pos)

    cell, _ = PointSymmetry.spglib_cell(lat, atoms, pos, mag_moments)
    # println(cell)
    ngrid, mapping, _grid = PointSymmetry.get_ir_reciprocal_mesh(cell, _kgrid_size, is_shift;
        is_time_reversal=is_time_reversal, symprec=PointSymmetry.SYMMETRY_TOLERANCE)

    _grid = Int.(_grid)
    grid = Vector{Vector{Int}}()
    for i in 1:length(mesh)
        k = _grid[3*(i-1)+1:3*i]
        # println(k)
        push!(grid, k)
    end
    # println(grid)

    # if size = [4, 4, 4]
    # grid = 
    # [[ 0  0  0]
    #  [ 1  0  0]
    #  [ 2  0  0]
    #  [-1  0  0]
    #  [ 0  1  0]
    #  [ 1  1  0]
    #  [ 2  1  0]
    #  [-1  1  0]
    #  ....      ]
    @assert grid[2][1] - grid[1][1] == 1 "expect the first index to iterate first"

    # mapping .+= 1
    kidx_unique = unique(mapping)

    # Convert irreducible k-points to DFTK conventions
    kirreds = [(is_shift ./ 2 .+ grid[ik]) ./ _kgrid_size for ik in kidx_unique]
    # println(kirreds)

    # Find the indices of the corresponding reducible k-points in `grid`, which one of the
    # irreducible k-points in `kirreds` generates.
    k_all_reducible = [findall(isequal(elem), mapping) for elem in kidx_unique]
    # println(k_all_reducible)

    # Number of reducible k-points represented by the irreducible k-point `kirreds[ik]`
    n_equivalent_k = length.(k_all_reducible)
    @assert sum(n_equivalent_k) == prod(kgrid_size)
    kweights = n_equivalent_k / sum(n_equivalent_k)
    # println("kweights: ", kweights)

    # This loop checks for reducible k-points, which could not be mapped to any irreducible
    # k-point yet even though spglib claims this can be done.
    # This happens because spglib actually fails for some non-ideal lattices, resulting
    # in *wrong results* being returned. See the discussion in
    # https://github.com/spglib/spglib/issues/101
    for (iks_reducible, k) in zip(k_all_reducible, kirreds), ikred in iks_reducible
        kred = (is_shift ./ 2 .+ grid[ikred]) ./ _kgrid_size

        found_mapping = any(symmetries) do symop
            # If the difference between kred and W' * k == W^{-1} * k
            # is only integer in fractional reciprocal-space coordinates, then
            # kred and S' * k are equivalent k-points
            all(isinteger, kred - (symop.S * k))
        end

        if !found_mapping
            error("The reducible k-point $kred could not be generated from " *
                  "the irreducible kpoints. This points to a bug in spglib.")
        end
    end

    new_map = zeros(Int, length(mesh))
    for (i, m) in enumerate(mapping)
        grid_address = grid[i]
        mapped_grid_address = grid[m]
        kidx = spglib_grid_address_to_index(mesh.mesh, grid_address)
        mapped_kidx = spglib_grid_address_to_index(mesh.mesh, mapped_grid_address)
        new_map[kidx] = mapped_kidx
    end

    # return kirreds, mapping, k_all_reducible
    return MeshMaps.MeshMap(new_map)

    # kcoords, kweights, symmetries = PointSymmetry.bzmesh_ir_wedge(_kgrid_size, symmetries; kshift=_kshift)
    # all_kcoords = PointSymmetry.unfold_kcoords(kcoords, symmetries)
    # Nk = reduce(*, kgrid_size)
    # @assert length(all_kcoords) == Nk
    # kindices = []
    # kmap = zeros(Int, Nk)
    # inv_kmap = Dict{Int,Vector{Int}}()
    # count = 0
    # for kpoint in kcoords
    #     all_kpoint = PointSymmetry.unfold_kcoords([kpoint,], symmetries)
    #     k0ind = _kcoords2ind(kpoint, kgrid_size, _kshift)
    #     push!(kindices, k0ind)
    #     inv_kmap[k0ind] = []
    #     count += length(all_kpoint)
    #     for k in all_kpoint
    #         kind = _kcoords2ind(k, kgrid_size, _kshift)
    #         kmap[kind] = k0ind
    #         push!(inv_kmap[k0ind], kind)
    #     end
    # end

    # @assert count == Nk

    # return MeshMaps.MeshMap(kindices, kmap, inv_kmap)
end








