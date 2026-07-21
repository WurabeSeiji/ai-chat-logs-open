import numpy as np

def simulate_steps(N, steps=10):
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
    u = np.random.randn(M)
    v = np.random.randn(M)
    u /= np.linalg.norm(u)
    v -= np.dot(u, v) * u
    v /= np.linalg.norm(v)
    
    R2 = 1.0
    s = 0.35
    X = np.sqrt(R2 + s**2) * u + 1j * s * v
    gamma = 0.1
    
    print(f"--- N={N} ---")
    for step in range(steps):
        theta = np.angle(X)
        S = np.sin(theta[None, :] - theta[:, None])
        K_tilde = A * S
        norm_K = np.linalg.norm(K_tilde, 2)
        if norm_K > 1e-12:
            K = K_tilde / norm_K
        else:
            K = K_tilde
            
        rank = np.linalg.matrix_rank(K)
        planes = rank // 2
        print(f"Step {step}: rank={rank}, planes={planes}")
        
        I = np.eye(M)
        U = np.linalg.inv(I - gamma * K) @ (I + gamma * K)
        X = U @ X

simulate_steps(4)
simulate_steps(5)
