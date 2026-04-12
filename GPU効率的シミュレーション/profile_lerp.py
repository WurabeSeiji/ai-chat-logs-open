import numpy as np
import time

def profile_lerp(N_side=100, S=10.0):
    N3 = N_side**3
    T = np.random.rand(11, 11, 11, 4).astype(np.float64)
    xi = np.random.rand(N3, 4) * 99.0
    
    t_start = time.time()
    
    # Block 1: coordinate math and floor
    t0 = time.time()
    p = xi[:, 1:4] / S
    idx0 = np.floor(p).astype(np.int32)
    idx0[:, 0] = np.clip(idx0[:, 0], 0, T.shape[0] - 2)
    idx0[:, 1] = np.clip(idx0[:, 1], 0, T.shape[1] - 2)
    idx0[:, 2] = np.clip(idx0[:, 2], 0, T.shape[2] - 2)
    idx1 = idx0 + 1
    t1 = time.time()
    
    # Block 2: Weights math
    w = p - idx0
    w0 = 1.0 - w
    w1 = w
    wx0 = w0[:, 0:1]; wy0 = w0[:, 1:2]; wz0 = w0[:, 2:3]
    wx1 = w1[:, 0:1]; wy1 = w1[:, 1:2]; wz1 = w1[:, 2:3]
    t2 = time.time()
    
    # Block 3: Memory fetch (8 corners)
    i0, j0, k0 = idx0[:, 0], idx0[:, 1], idx0[:, 2]
    i1, j1, k1 = idx1[:, 0], idx1[:, 1], idx1[:, 2]
    
    c000 = T[i0, j0, k0]
    c100 = T[i1, j0, k0]
    c010 = T[i0, j1, k0]
    c110 = T[i1, j1, k0]
    c001 = T[i0, j0, k1]
    c101 = T[i1, j0, k1]
    c011 = T[i0, j1, k1]
    c111 = T[i1, j1, k1]
    t3 = time.time()
    
    # Block 4: Weight multiplication
    w000 = wx0 * wy0 * wz0
    w100 = wx1 * wy0 * wz0
    w010 = wx0 * wy1 * wz0
    w110 = wx1 * wy1 * wz0
    w001 = wx0 * wy0 * wz1
    w101 = wx1 * wy0 * wz1
    w011 = wx0 * wy1 * wz1
    w111 = wx1 * wy1 * wz1
    t4 = time.time()
    
    # Block 5: Final sum
    X_interp = (c000 * w000 + c100 * w100 + 
                c010 * w010 + c110 * w110 + 
                c001 * w001 + c101 * w101 + 
                c011 * w011 + c111 * w111)
    t5 = time.time()
    
    print(f"Total time   : {t5 - t_start:.4f} s")
    print(f"1. Floor/Clip: {t1 - t0:.4f} s")
    print(f"2. Weight    : {t2 - t1:.4f} s")
    print(f"3. ArrayFetch: {t3 - t2:.4f} s")
    print(f"4. WeightMult: {t4 - t3:.4f} s")
    print(f"5. Final Sum : {t5 - t4:.4f} s")

if __name__ == '__main__':
    profile_lerp()
