import numpy as np

for N in range(3, 10):
    M = N * (N - 1) // 2
    B = np.zeros((N, M))
    col = 0
    for i in range(N):
        for j in range(i+1, N):
            B[i, col] = 1
            B[j, col] = 1
            col += 1
            
    np.random.seed(42 + N)
    theta = np.random.rand(M) * 2 * np.pi
    u = np.sin(theta)
    v = np.cos(theta)
    
    W = np.zeros((M, 2*N))
    for k in range(N):
        W[:, 2*k] = B[k, :] * u
        W[:, 2*k+1] = B[k, :] * v
        
    rank_W = np.linalg.matrix_rank(W)
    
    A = B.T @ B - 2 * np.eye(M)
    S = np.sin(theta[None, :] - theta[:, None])
    K = A * S
    rank_K = np.linalg.matrix_rank(K)
    
    print(f"N={N}, M={M}, rank(W)={rank_W}, rank(K)={rank_K}")
