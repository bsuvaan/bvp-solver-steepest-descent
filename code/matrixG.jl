### A Pluto.jl notebook ###
# v0.19.47

using Markdown
using InteractiveUtils

# ╔═╡ 28fd47c6-cf52-482d-96fe-019007e7bb30
using LinearAlgebra

# ╔═╡ 2969a8ea-033f-45af-beaf-79be6a0e1c5a
function as_square_matrix(vv::Vector{<:AbstractVector})
	    M = length(vv)
	    @assert all(length(v) == M for v in vv) "длины под-векторов должны равняться $M"
		vv_matrix = reshape(vcat(vv...), M, M) 
	
	    return   vv_matrix'# столбцы = исходные векторы
end

# ╔═╡ d657014c-317f-11f0-1204-6dbe52bf6751
function p(x::Float64)
	return 1.0 + x
end

# ╔═╡ 65dea137-fd6f-4c8f-af6c-4a7d6906ae22
function p_der(x::Float64)
	return 1
end

# ╔═╡ 3a76ac7e-a012-4e75-bbb6-e3328862c2e3
function q(x::Float64)
	return 1
end

# ╔═╡ 6b34b397-d97b-41d7-bd07-ec1346d8a126
function u(x::Float64)
	return x*(1-x)+(sin(pi*x))^2
end

# ╔═╡ 2f2c599c-bd92-4151-b4d1-3564e3acbc22
function u_der(x::Float64)
	return 1.0 - 2.0*x + pi*sin(2*pi*x)
end

# ╔═╡ 3f473c08-5bf5-4835-8e01-66ac1147e835
function u_der_der(x::Float64)
	return 2*pi^2 * cos(2*pi*x) - 2
end

# ╔═╡ a7dd712e-f8dc-4b8d-8664-269635127270
function f(x::Float64)
	f = p_der(x) * u_der(x) + p(x) * u_der_der(x) - q(x)*u(x)
	return f
end

# ╔═╡ a522f94d-dca4-4b9a-b7a6-ce7b84e11829
# Здесь я задаю параметры сетки 
begin
	x0 = 0
	x_M = 1
	M = 100
	Δx = abs(x_M-x0)/M
	#X = range(x0, stop = x_M, length = M)

	u_0 = 0
	u_1 = 0
end

# ╔═╡ 7b62f5d8-532a-417c-b26a-349b8a4b06c4
function P_plus(i::Int)
	return p(i* Δx + Δx/2)
end

# ╔═╡ 79398937-e2f2-48a1-bc0c-d297bfe4ef77
function P_minus(i)
	return p(i* Δx - Δx/2)
end

# ╔═╡ 35b32ce7-48d8-4e99-8ef4-8d99de824a83
begin
	U_analyt = u.(x0:Δx:x_M)
end

# ╔═╡ dc503cdf-13b4-4d02-a951-50665cb003d0
# строим вектор F
begin
	F = []
	for i in 1:(M-1)
		f_i = - Δx^2 * f(Δx*i)
		push!(F, f_i)
	end
end

# ╔═╡ aeec82cc-8efa-4601-b980-84dd437052ab
# строим вектор G
begin
	G_vec = Vector{Vector{Float64}}()
	for i in 1:(M-1)
		G_row = zeros(M-1)::Vector{Float64}
		#k = i - 1 # курсор 
		if i == 1
			#continue
			#G_row[1] = 1.0
			G_row[i] = P_minus(i) + P_plus(i) + q(i*Δx) * Δx^2
			G_row[i+1] = -1.0 * P_plus(i)
		elseif i == (M-1)
			#continue
			#G_row[M] = 1.0
			G_row[i-1] = -1.0 * P_minus(i)
			G_row[i] = P_minus(i) + P_plus(i) + q(i*Δx) * Δx^2
		else
			x = i * Δx
			G_row[i-1] = -1.0 * P_minus(i)
			G_row[i] = P_minus(i) + P_plus(i) + q(i*Δx) * Δx^2
			G_row[i+1] = -1.0 * P_plus(i)
		end
		push!(G_vec, G_row)
	end	
end

# ╔═╡ 97a044f2-91d5-4b93-82b6-d159cc87ccca
length(G_vec)

# ╔═╡ 9b4dbee6-ee97-49b7-b898-1c61fe7f166e
G = as_square_matrix(G_vec)

# ╔═╡ 45b32a13-16ae-4a9d-a25f-d1e45d699a61
function find_nonsymmetric_indices(A::LinearAlgebra.Adjoint{Float64, Matrix{Float64}})
    n, m = size(A)
    @assert n == m "Матрица должна быть квадратной"
    mismatches = []
    for i in 1:n, j in i+1:n
        if A[i, j] != A[j, i]
            push!(mismatches, A[i, j] - A[j, i])
        end
    end
    return mismatches
end

# ╔═╡ 149f4724-aa0d-4b59-bb30-794c789d3470
find_nonsymmetric_indices(G)

# ╔═╡ d27bb7b5-5517-4ba1-b046-3f0b5fc59c71
G

# ╔═╡ c60b8a95-535c-464c-9a3c-41c290582826
typeof(F)

# ╔═╡ 4eb00e5d-e090-45b9-a5ce-4ab496e40fcd
λmin = (pi^2 * p(0.0) + q(0.0)) * Δx^2

# ╔═╡ 5302aa93-d7c2-4a4c-b91a-941e635c6a59
λmax = 4*p(Δx) + q(Δx)*Δx^2

# ╔═╡ 49048de6-e94f-4e9b-9601-47e24b1af66e
function gerschgorin_bounds(A::LinearAlgebra.Adjoint{Float64, Matrix{Float64}})
    n, m = size(A)
    @assert n == m "Матрица должна быть квадратной"
    centers = diag(A)
    radii = similar(centers)
    lowers = similar(centers)
    uppers = similar(centers)

    for i in 1:n
        # сумма модулей внедиагональных элементов строки i
        radii[i] = sum(abs.(A[i, j] for j in 1:n if j != i))
        lowers[i] = centers[i] - radii[i]
        uppers[i] = centers[i] + radii[i]
    end

    λ_min = minimum(lowers)
    λ_max = maximum(uppers)

    # Соберем результат
    bounds = [(centers[i], radii[i], lowers[i], uppers[i]) for i in 1:n]
	alpha = (2.0) / (λmin + λmax)
    return alpha
end



# ╔═╡ 07caa872-caa8-4372-91d6-2bcc18ec5601
gerschgorin_bounds(G)

# ╔═╡ c39cf8c7-166a-4698-871e-9b8f11d52ad1
function optimum_step()
	λmin = (pi^2 * p(0.0) + q(0.0)) * Δx^2
	λmax = 4*p(Δx) + q(Δx)*Δx^2
	return alpha = (2.0) / (λmin + λmax)
end

# ╔═╡ d0a5dfe3-8bf5-4a39-adf8-26c6ef8d8880
optimum_step()

# ╔═╡ ab7e18f8-1fa3-4acb-9af7-88c24aad848e
function steepest_descent(G, F;
                        y0=zeros(length(F)),
                        tol=1e-6, maxiter=50_000)

    #α = optimum_step()         # оптимальный постоянный шаг
    y = copy(y0)
    r = F - G*y #невязка
    res2 = dot(r, r) #квадрат невязки 
	res_hist = []
    hist = [y]

    for k in 1:maxiter
		res = sqrt(res2) # норма невязки 
		if res <= tol
			return hist, res_hist
		else
			r = F-G*y              # пересчитываем невязку
        	res2 = dot(r,r)
			z = G*r
			α  = res2 / dot(r, z) 
			y = y + α .* r             # y_{k+1} = y_k + α * r_k

			push!(res_hist, sqrt(res2))
        	push!(hist, y)
			#print(y, "Невязка ", res)
		end
    end
    # @warn "Достигнут maxiter; ‖r‖∞ = $res"res
    return hist, res_hist
end


# ╔═╡ d08a79fa-45b4-4cd9-a98c-51aabaf67c83
begin 
	D = [1.67 0.32  0.12 0.57;
		 0.32 4.17 0.65 0.15;
		 0.12 0.65 3.15 0.22;
		 0.57 0.15 0.22 1.84]
	b = [1.34;
		0.85; 
		1.29; 
		2.11]
	eps = 0.001
	steepest_descent(D, b;
                        y0=zeros(length(b)),
                        tol=eps)
end

# ╔═╡ a4051e5c-1a9b-4b1d-a0e2-9ef1b77349e1
F

# ╔═╡ 53014726-fcfc-4647-a7b1-e8f103d4518c
begin
	Y = steepest_descent(G, F)
end

# ╔═╡ 217c5e0e-1595-405c-a82f-e529868c2f7c
U_analyt

# ╔═╡ e11a6382-24ce-42fe-ab3e-c4dd7185bdbd
norm([0.0
1.21415e17
1.094e16
-5.76745e17
-3.94915e17
1.16201e18
1.08537e18
-1.04091e18
-7.47217e17
0.0
], 2) 

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.0"
manifest_format = "2.0"
project_hash = "ac1187e548c6ab173ac57d4e72da1620216bce54"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"
"""

# ╔═╡ Cell order:
# ╠═28fd47c6-cf52-482d-96fe-019007e7bb30
# ╠═2969a8ea-033f-45af-beaf-79be6a0e1c5a
# ╠═d657014c-317f-11f0-1204-6dbe52bf6751
# ╠═65dea137-fd6f-4c8f-af6c-4a7d6906ae22
# ╠═7b62f5d8-532a-417c-b26a-349b8a4b06c4
# ╠═79398937-e2f2-48a1-bc0c-d297bfe4ef77
# ╠═3a76ac7e-a012-4e75-bbb6-e3328862c2e3
# ╠═6b34b397-d97b-41d7-bd07-ec1346d8a126
# ╠═2f2c599c-bd92-4151-b4d1-3564e3acbc22
# ╠═3f473c08-5bf5-4835-8e01-66ac1147e835
# ╠═a7dd712e-f8dc-4b8d-8664-269635127270
# ╠═a522f94d-dca4-4b9a-b7a6-ce7b84e11829
# ╠═35b32ce7-48d8-4e99-8ef4-8d99de824a83
# ╠═dc503cdf-13b4-4d02-a951-50665cb003d0
# ╠═aeec82cc-8efa-4601-b980-84dd437052ab
# ╠═97a044f2-91d5-4b93-82b6-d159cc87ccca
# ╠═9b4dbee6-ee97-49b7-b898-1c61fe7f166e
# ╠═45b32a13-16ae-4a9d-a25f-d1e45d699a61
# ╠═149f4724-aa0d-4b59-bb30-794c789d3470
# ╠═d27bb7b5-5517-4ba1-b046-3f0b5fc59c71
# ╠═c60b8a95-535c-464c-9a3c-41c290582826
# ╠═49048de6-e94f-4e9b-9601-47e24b1af66e
# ╠═07caa872-caa8-4372-91d6-2bcc18ec5601
# ╠═4eb00e5d-e090-45b9-a5ce-4ab496e40fcd
# ╠═5302aa93-d7c2-4a4c-b91a-941e635c6a59
# ╠═c39cf8c7-166a-4698-871e-9b8f11d52ad1
# ╠═d0a5dfe3-8bf5-4a39-adf8-26c6ef8d8880
# ╠═ab7e18f8-1fa3-4acb-9af7-88c24aad848e
# ╠═d08a79fa-45b4-4cd9-a98c-51aabaf67c83
# ╠═a4051e5c-1a9b-4b1d-a0e2-9ef1b77349e1
# ╠═53014726-fcfc-4647-a7b1-e8f103d4518c
# ╠═217c5e0e-1595-405c-a82f-e529868c2f7c
# ╠═e11a6382-24ce-42fe-ab3e-c4dd7185bdbd
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
