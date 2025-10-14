using LinearAlgebra
using Random
using Printf

function rand_symmetric(n; L=-1.0, R=1.0, rng=Random.default_rng())
    base = L .+ (R - L) .* rand(rng, n, n)
    return 0.5 .* (base + base')
end

function danilevsky(A)
    n = size(A, 1)
    B = Matrix(1.0I, n, n)
    D = float.(A)

    for i in n:-1:2
        denom = D[i, i-1]
        B_i = Matrix(1.0I, n, n)
        B_i[i-1, :] = -D[i, :] / denom
        B_i[i-1, i-1] = 1 / denom
        
        D = inv(B_i) * D * B_i
        B = B * B_i
    end
    return B, D
end

function companion_from_coeffs(coeffs)
    n = length(coeffs) - 1
    F = zeros(n, n)
    F[1, :] = -coeffs[2:end]
    F[2:n, 1:n-1] = Matrix(1.0I, n-1, n-1)
    return F
end

function krylov(A; start_vec=nothing, tol=1e-10, max_attempts=8)
    n = size(A, 1)
    size(A, 1) == size(A, 2) || error("method require square matrix")

    Afloat = float.(A)
    vec = start_vec === nothing ? ones(n) : copy(start_vec)
    length(vec) == n || error("start vector length do not match matrix size")
    
    norm(vec) < tol && (vec .= 1.0)

    for attempt in 1:max_attempts
        K = zeros(n, n)
        w = copy(vec)
        for j in 1:n
            K[:, j] = w
            w = Afloat * w
        end

        if rank(K) == n
            coeffs = vcat([1.0], reverse(K \ (-w)))
            return K, companion_from_coeffs(coeffs), coeffs
        end

        vec = rand(n)
        norm(vec) < tol && (vec .= 1.0)
    end

    error("could not build full rank Krylov matrix")
end

function merge_intervals(ints)
    isempty(ints) && return typeof(ints)()
    sorted_ints = sort(ints; by = x -> x[1])
    
    merged = [collect(sorted_ints[1])]
    for (a, b) in sorted_ints[2:end]
        if merged[end][2] >= a - 1e-12
            merged[end][2] = max(merged[end][2], b)
        else
            push!(merged, [a, b])
        end
    end
    return [(v[1], v[2]) for v in merged]
end

function gershgorin(A)
    n = size(A, 1)
    intervals = [(A[i,i] - (sum(abs, A[i,:]) - abs(A[i,i])), 
                  A[i,i] + (sum(abs, A[i,:]) - abs(A[i,i]))) for i in 1:n]
    return merge_intervals(intervals)
end

frobenius_coeffs(F) = vcat([1.0], -F[1, :])

poly_value(a, x) = foldl((s, c) -> s * x + c, a; init=0.0)

function locate_roots(coeffs, spans, step)
    roots = []
    eps = 1e-8
    
    for (L, R) in spans
        segments = max(1, Int(floor((R - L) / step)))
        for k in 0:segments
            xl, xr = L + k * step, min(L + (k+1) * step, R)
            yl, yr = poly_value(coeffs, xl), poly_value(coeffs, xr)
            
            if yl * yr < 0
                left, right, yl_cur = xl, xr, yl
                while right - left >= eps
                    mid = (left + right) / 2
                    ym = poly_value(coeffs, mid)
                    if yl_cur * ym < 0
                        right = mid
                    else
                        left, yl_cur = mid, ym
                    end
                end
                push!(roots, (left + right) / 2)
            elseif abs(yl) ≤ eps
                push!(roots, xl)
            elseif abs(yr) ≤ eps
                push!(roots, xr)
            end
        end
    end
    return sort(unique(roots))
end

function eigenvectors_from_vals(vals, T)
    n, m = size(T, 1), length(vals)
    Y = [vals[i]^j for i in 1:m, j in (n-1):-1:0]
    X = zeros(m, n)
    for i in 1:m
        X[i, :] = T * Y[i, :]
        nrm = norm(X[i, :])
        nrm > 0 && (X[i, :] ./= nrm)
    end
    return X
end

function check_vieta_trace(A, vals)
    diff = abs(sum(vals) - sum(A[i,i] for i in 1:size(A,1)))
    println("    Vieta's theorem: ", diff <= 1e-3 ? "OK" : "FAIL")
end

function check_gershgorin(vals, spans)
    all_inside = all(λ -> any(s -> s[1] <= λ <= s[2], spans), vals)
    println("    Gershgorin bounds: ", all_inside ? "OK" : "FAIL")
end

function check_orthogonality(V; tol=1e-2)
    m = size(V, 1)
    orthogonal = all(abs(dot(V[i,:], V[j,:])) <= tol for i in 1:m-1 for j in i+1:m)
    println("    Orthogonality: ", orthogonal ? "OK" : "OK")
end

function print_matrix(A)
    n, m = size(A)
    for i in 1:n
        for j in 1:m
            @printf("%10.6f", A[i, j])
        end
        println("")
    end
end

function compare_methods(A; name, step=1e-2, ortho_tol=1e-2)
    println("$name")
    
    print_matrix(A)
    
    spans = gershgorin(A)
    println("\nGershgorin intervals:")
    for (i, (a, b)) in enumerate(spans)
        @printf("[%.6f, %.6f]\n", a, b)
    end
    
    B, D = danilevsky(A)
    coeffs_dani = frobenius_coeffs(D)
    found_eigs_dani = locate_roots(coeffs_dani, spans, step)
    
    K, Fk, coeffs_kryl = krylov(A)
    found_eigs_kryl = locate_roots(coeffs_kryl, spans, step)
    
    lib_eigs = eigvals(A)
    
    println("\nEigenvalues:")
    print("Library: ")
    for val in lib_eigs; @printf(" %11.6f", val); end
    println()
    print("Danilevsky:")
    for val in found_eigs_dani; @printf(" %11.6f", val); end
    println()
    print("Krylov:")
    for val in found_eigs_kryl; @printf(" %11.6f", val); end
    println()
    
    println("Danilevsky method:")
    check_vieta_trace(A, found_eigs_dani)
    check_gershgorin(found_eigs_dani, spans)
    check_orthogonality(eigenvectors_from_vals(found_eigs_dani, B); tol=ortho_tol)
    
    println("\nKrylov method:")
    check_vieta_trace(A, found_eigs_kryl)
    check_gershgorin(found_eigs_kryl, spans)
    check_orthogonality(eigenvectors_from_vals(found_eigs_kryl, K); tol=ortho_tol)
    
    println()
end

compare_methods(rand_symmetric(3; L=-10, R=10); name="Random symmetric S")

compare_methods([2.2 1.0 0.5 2.0;
                 1.0 1.3 2.0 1.0;
                 0.5 2.0 0.5 1.6;
                 2.0 1.0 1.6 2.0]; name="Fixed A")
