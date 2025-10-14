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

danilevsky() = begin
    A = [1.0 2.0 3.0; 2.0 4.0 5.0; 3.0 5.0 6.0]
    @test is_symmetric(A)

    n = size(A, 1)
    b = [1.0, 2.0, 3.0]

    @test A[n, n-1] ≠ 0

    D = copy(A)
    for row in n:-1:2
        if D[row, row-1] == 0
            # TODO: check if it is possible to swap rows and columns 
        end
        B = get_B(D, row)

        C = D * B
        D = B \ C
    end
    println(D)
    return D
end

danilevsky()
