import numpy as np

def simulate(N):
    M = N * (N - 1) // 2
    B = np.zeros((N, M))
    col = 0
    for i in range(N):
        for j in range(i+1, N):
            B[i, col] = 1
            B[j, col] = 1
            col += 1
            
    A = B.T @ B - 2 * np.eye(M)
    
    np.random.seed(42 + N)
    # Initialize state as in the paper
    u = np.random.randn(M)
    v = np.random.randn(M)
    u /= np.linalg.norm(u)
    v -= np.dot(u, v) * u
    v /= np.linalg.norm(v)
    
    R2 = 1.0
    s = 0.35
    X0 = np.sqrt(R2 + s**2) * u + 1j * s * v
    
    theta = np.angle(X0)
    
    S = np.sin(theta[None, :] - theta[:, None])
    K_tilde = A * S
    
    norm_K = np.linalg.norm(K_tilde, 2)
    if norm_K > 1e-12:
        K = K_tilde / norm_K
    else:
        K = K_tilde
        
    rank = np.linalg.matrix_rank(K)
    planes = rank // 2
    print(f"N={N}, M={M}, Generator Rank={rank}, Rotational Planes={planes}")

for n in range(3, 8):
    simulate(n)
