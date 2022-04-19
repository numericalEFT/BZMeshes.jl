@testset "Base Mesh" begin

    function dispersion(k)
        me = 0.5
        return dot(k, k) / 2me
    end

	  function density(k)
        T = 0.01
        μ = 1.0

        ϵ = dispersion(k) - μ

        # return 1 / (exp((ϵ) / T) + 1.0)
        return (π*T)^2/((π*T)^2 + ϵ^2)
        # return (exp((ϵ) / T) / T)/(exp((ϵ) / T) + 1.0)^2
    end

    @testset "Uniform Mesh" begin
        δ = 1e-3

        # index convert
        N, DIM = 4, 2
        for i in 1:N^(DIM)
            # println("$i -> $(BaseMesh._ind2inds(i, N, DIM))")
            @test i == BaseMesh._inds2ind(BaseMesh._ind2inds(i, N, DIM), N)
        end

        # 2D
        N, DIM = 4, 2
        origin = [0.0 0.0]
        latvec = [2 0;1 sqrt(3)]'
        # latvec = [1 0; 0 1]'

        umesh = UniformMesh{DIM, N}(origin, latvec)

        for i in 1:length(umesh)
            for j in 1:DIM
                shift = zeros(Float64, DIM)
                indshift = zeros(Int, DIM)

                inds = BaseMesh._ind2inds(i, N, DIM)
                # println(inds)

                shift = δ .* latvec[:, j]
                indshift[j] = 0
                for k in 1:DIM
                    if inds[k] == N
                        inds[k] = N-1
                        if k == j
                            indshift[j] = 0
                        end
                    end
                end
                ind = BaseMesh._inds2ind( inds + indshift, N)
                # @test ind == floor(umesh, umesh[i] + shift)
                @test BaseMesh._ind2inds(ind, N, DIM) == BaseMesh._ind2inds(floor(umesh, umesh[i] + shift), N, DIM)

                inds = BaseMesh._ind2inds(i, N, DIM)
                shift = -δ .* latvec[:, j]
                indshift[j] = -1
                for k in 1:DIM
                    if inds[k] == N
                        inds[k] = N-1
                        if k == j
                            indshift[j] = 0
                        end
                    end
                end
                if inds[j] == 1
                    indshift[j] = 0
                end
                ind = BaseMesh._inds2ind( inds + indshift, N)
                # @test ind == floor(umesh, umesh[i] + shift)
                @test BaseMesh._ind2inds(ind, N, DIM) == BaseMesh._ind2inds(floor(umesh, umesh[i] + shift), N, DIM)
            end
        end

        # 3D
        N, DIM = 3, 3
        origin = [0.0 0.0 0.0]
        latvec = [1.0 0 0; 0 1.0 0; 0 0 1.0]'

        umesh = UniformMesh{DIM, N}(origin, latvec)

        for i in 1:length(umesh)
            for j in 1:DIM
                shift = zeros(Float64, DIM)
                indshift = zeros(Int, DIM)

                inds = BaseMesh._ind2inds(i, N, DIM)
                # println(inds)

                shift = δ .* latvec[:, j]
                indshift[j] = 0
                for k in 1:DIM
                    if inds[k] == N
                        inds[k] = N-1
                        if k == j
                            indshift[j] = 0
                        end
                    end
                end
                ind = BaseMesh._inds2ind( inds + indshift, N)
                # @test ind == floor(umesh, umesh[i] + shift)
                @test BaseMesh._ind2inds(ind, N, DIM) == BaseMesh._ind2inds(floor(umesh, umesh[i] + shift), N, DIM)

                inds = BaseMesh._ind2inds(i, N, DIM)
                shift = -δ .* latvec[:, j]
                indshift[j] = -1
                for k in 1:DIM
                    if inds[k] == N
                        inds[k] = N-1
                        if k == j
                            indshift[j] = 0
                        end
                    end
                end
                if inds[j] == 1
                    indshift[j] = 0
                end
                ind = BaseMesh._inds2ind( inds + indshift, N)
                # @test ind == floor(umesh, umesh[i] + shift)
                @test BaseMesh._ind2inds(ind, N, DIM) == BaseMesh._ind2inds(floor(umesh, umesh[i] + shift), N, DIM)
            end
        end

    end

    @testset "Interp and Integral" begin
        N, DIM = 200, 2
        origin = [0.0 0.0]
        latvec = [π 0; 0 π]'
        umesh = UniformMesh{DIM, N}(origin, latvec)

        # f(x) = x[1] + 2 * x[2] + x[1] * x[2]
        f(x) = sin(x[1]) + cos(x[2])

        data = zeros(Float64, size(umesh))

        for i in 1:length(umesh)
            data[i] = f(umesh[i])
        end

        # interpolate
        testN = 10
        xlist = rand(testN) * π
        ylist = rand(testN) * π
        for x in xlist
            for y in ylist
                @test isapprox(f([x, y]), interp(data, umesh, [x, y]), rtol = 2e-2)
            end
        end

        # integrate
        integral = integrate(data, umesh)
        # println("integral=$(integral)")
        @test isapprox(integral, 2π, rtol = 1e-3)
    end

end
