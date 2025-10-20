using Test
using Random
using LinearAlgebra

Random.seed!(42)

generate_matrix(n) = begin
    A = rand(n, n)
    for i in 1:n
        A[i, i] = sum(abs(A[i, j]) for j in 1:n if j != i) + rand() + 1.0
    end
    return A
end

check_convergence(A) = begin
    for i in axes(A, 1)
        if abs(A[i, i]) <= sum(abs(A[i, j]) for j in axes(A, 2) if j != i)
            return false
        end
    end
    return true
end

yakobi(A, b) = begin
    if !check_convergence(A)
        error("Matrix is not convergent")
    end

    x = zeros(size(A, 1))
    eps = 1e-6
    max_iterations = 10000
    for i in 1:max_iterations
        x_prev = copy(x)
        for j in axes(A, 1)
            x[j] = (b[j] - sum(A[j, k] * x_prev[k] for k in axes(A, 2) if k != j)) / A[j, j]
        end
        if norm(x - x_prev) < eps
            println("Desired accuracy achieved")
            return x, i
        end
    end
    println("Desired accuracy not achieved")
    return x, max_iterations
end

zeidel(A, b) = begin
    if !check_convergence(A)
        error("Matrix is not convergent")
    end

    x = zeros(size(A, 1))
    eps = 1e-6
    max_iterations = 10000
    for i in 1:max_iterations
        x_prev = copy(x)
        for j in axes(A, 1)
            x[j] = (b[j] - sum(A[j, k] * x[k] for k in axes(A, 2) if k != j)) / A[j, j]
        end
        if norm(x - x_prev) < eps
            println("Desired accuracy achieved")
            return x, i
        end
    end
    println("Desired accuracy not achieved")
    return x, max_iterations
end

@test zeidel([4.0 0.24 -0.08; 0.09 3.0 -0.15; 0.04 -0.08 4.0], [8.0, 9.0, 20.0])[1] ≈ [1.909, 3.195, 5.045] atol=1e-3
@test yakobi([4.0 0.24 -0.08; 0.09 3.0 -0.15; 0.04 -0.08 4.0], [8.0, 9.0, 20.0])[1] ≈ [1.909, 3.195, 5.045] atol=1e-3

main() = begin
    A = generate_matrix(500)
    b = rand(500)
    
    time_zeidel = @elapsed begin
        x, iterations = zeidel(A, b)
    end
    println("Zeidel:")
    println("Iterations: ", iterations)
    println("Time: ", time_zeidel * 1000, " ms")
    println("Error: ", norm(A * x - b))
    
    time_yakobi = @elapsed begin
        x, iterations = yakobi(A, b)
    end
    println("Yakobi:")
    println("Iterations: ", iterations)
    println("Time: ", time_yakobi * 1000, " ms")
    println("Error: ", norm(A * x - b))
end

main()