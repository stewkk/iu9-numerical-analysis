using Test
using Random
using LinearAlgebra

Random.seed!(42)

generate_symmetric_matrix(n) = begin
    A = rand(n, n)
    A = A + A'
    return A
end

is_symmetric(A) = begin
    return A ≈ A'
end

@test is_symmetric(generate_symmetric_matrix(5))

get_B(A, i) = begin
    n = size(A, 1)
    B = Matrix{Float64}(I, n, n)
    for j in 1:n
        if j != i-1
            B[i-1, j] = - A[i, j] / A[i,i-1]
        end
    end
    B[i-1, i-1] = 1/A[i, i-1]
    return B
end

get_frobenius_matrix(A) = begin
    n = size(A, 1)
    D = copy(A)
    B_res = Matrix{Float64}(I, n, n)
    for row in n:-1:2
        if D[row, row-1] == 0
            println("swapping")
            m = row-2
            while m > 0 && D[row, m] == 0
                m -= 1
            end
            if m == 0
                error("")
            end
            D[:, [m, row-1]] = D[:, [row-1, m]]
            D[m, :], D[row-1, :] = D[row-1, :], D[m, :]
        end
        B = get_B(D, row)

        C = D * B
        D = B \ C
        B_res = B_res * B
    end
    return D, B_res
end

@test get_frobenius_matrix([2.2 1.0 0.5 2.0; 1.0 1.3 2.0 1.0; 0.5 2.0 0.5 1.6; 2.0 1.0 1.6 2.0]) == ([6.0 0.20000000000000182 -12.735000000000005 2.7616; 1.0 0.0 0.0 0.0; -1.6037505170817941e-18 1.0 1.1456124527020951e-17 -8.040135925636728e-18; 2.115006742527263e-17 -1.121575634346075e-16 0.9999999999999999 1.4384720710653878e-16], [-0.23112480739599386 1.078582434514638 1.6510015408320498 -1.1587057010785826; 0.081243871690713 -0.1367138254657515 -1.6409581173833871 -0.2739095111360136; 0.2381285894382967 -1.2627819022272029 -0.41315310267544514 0.3695755708082368; 0.0 0.0 0.0 1.0])

find_gershgorin_intervals(A) = begin
    n = size(A, 1)
    intervals = []
    for i in 1:n
        interval = (A[i, i] - sum(abs(A[i, j]) for j in 1:n if j != i), A[i, i] + sum(abs(A[i, j]) for j in 1:n if j != i))
        push!(intervals, interval)
    end
    return intervals
end

@test find_gershgorin_intervals([-2.0 0.5 0.5; -0.5 -3.5 1.5; 0.8 -0.5 0.5]) == [(-3.0, -1.0), (-5.5, -1.5), (-0.8, 1.8)]

find_polynom_roots(coeffs, intervals) = begin
    eval_poly(x) = begin
        result = coeffs[1]
        for i in 2:lastindex(coeffs)
            result = result * x + coeffs[i]
        end
        return result
    end
    
    sorted_intervals = sort(intervals, by = x -> x[1])
    merged = [sorted_intervals[1]]
    for (a, b) in sorted_intervals[2:end]
        last_a, last_b = merged[end]
        if a <= last_b
            merged[end] = (last_a, max(b, last_b))
        else
            push!(merged, (a, b))
        end
    end
    
    roots = []
    tol = 1e-5
    max_iter = 1000
    intervals_count = 100
    
    for (interval_start, interval_end) in merged
        step = (interval_end - interval_start) / intervals_count
        
        for i in 0:(intervals_count-1)
            left = interval_start + i * step
            right = interval_start + (i + 1) * step
            f_left = eval_poly(left)
            f_right = eval_poly(right)
            
            if f_left * f_right < 0
                a, b = left, right
                fa, fb = f_left, f_right
                
                for iter in 1:max_iter
                    mid = (a + b) / 2
                    f_mid = eval_poly(mid)
                    
                    if abs(f_mid) < tol || (b - a) / 2 < tol
                        push!(roots, mid)
                        break
                    end
                    
                    if fa * f_mid < 0
                        b = mid
                        fb = f_mid
                    else
                        a = mid
                        fa = f_mid
                    end
                end
            end
        end
        
    end
    
    unique_roots = []
    for root in sort(roots)
        if isempty(unique_roots) || abs(root - unique_roots[end]) > tol
            push!(unique_roots, round(root, digits=4))
        end
    end
    
    return unique_roots
end

@test find_polynom_roots([1.0, -6.0, -0.2, 12.735, -2.7616], [(-3.0, -1.0), (0.0, 1.0), (1.0, 2.0), (4.0, 6.0)]) == [-1.4201, 0.2226, 1.5454, 5.652]

danilevsky(A) = begin
    @test is_symmetric(A)
    n = size(A, 1)

    P, B = get_frobenius_matrix(A)
    first_row = P[1, :]
    first_row *= -1
    first_row = vcat(1, first_row)
    
    gershgorin_intervals = find_gershgorin_intervals(A)

    eigenvalues = find_polynom_roots(first_row, gershgorin_intervals)
    
    eigenvectors = []
    for eigenvalue in eigenvalues
        eigenvector = B * [eigenvalue^i for i in n-1:-1:0]
        push!(eigenvectors, eigenvector)
    end

    return eigenvalues, eigenvectors
end

@test danilevsky([1.0 2.0 3.0; 2.0 4.0 5.0; 3.0 5.0 6.0]) == ([-0.5157, 0.1709, 11.3448], [[-1.246840957692308, -0.5550354253846155, 1.0], [1.8018743038461542, -2.2469445823076923, 1.0], [0.4451049846153863, 0.8018970092307711, 1.0]])
A = [1.0 2.0 3.0; 2.0 3.0 0.0; 0.0 0.0 6.0]
@test danilevsky(A + A') == ([-0.9841, 7.9557, 13.0284], [[2.4788210675, -4.328033333333333, 1.0], [-2.757219792500001, -1.3481, 1.0], [0.19513387999999776, 0.34279999999999955, 1.0]])




