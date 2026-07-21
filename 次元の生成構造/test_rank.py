import numpy as np

for N in range(3, 10):
    M = N * (N - 1) // 2
    # Create incidence matrix B
    B = np.zeros((N, M))
    col = 0
    for i in range(N):
        for j in range(i+1, N):
            B[i, col] = 1
            B[j, col] = 1
            col += 1
            
    A = B.T @ B - 2 * np.eye(M)
    
    # Random phases
    np.random.seed(42 + N)
    theta = np.random.rand(M) * 2 * np.pi
    
    S = np.sin(theta[None, :] - theta[:, None])
    K = A * S
    
    # rank of K
    rank = np.linalg.matrix_rank(K)
    print(f"N={N}, M={M}, rank(K)={rank}")
