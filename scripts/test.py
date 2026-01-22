import numpy as np
import matplotlib.pyplot as plt
from numba import njit

@njit
def set_value(A, p1, p2, X, L):
    for i in range(L):
        for j in range(L):
            x = i + 0.5
            y = j + 0.5
            rx = np.abs(x - X[p1,0])
            rx = np.minimum(rx, L-rx)
            ry = np.abs(y - X[p1,1])
            ry = np.minimum(ry, L-ry)

            delta_1 = kernel(rx)*kernel(ry)
            rx = np.abs(x - X[p2,0])
            rx = np.minimum(rx, L-rx)
            ry = np.abs(y - X[p2,1])
            ry = np.minimum(ry, L-ry)

            delta_2 = kernel(rx)*kernel(ry)
            A[p1, p2] += delta_1 * delta_2

@njit
def kernel(x):
    # if x < 1:
        # return 1-x
    # else:
        # return 0
    if x <= 0.5:
        return 1/3*(1 + np.sqrt(1 - 3*x*x))
    elif x <= 1.5:
        return 1/6*(5 - 3*x - np.sqrt(-2 + 6*x - 3*x*x))
    else:
        return 0

if __name__ == "__main__":

    N_repeat = 1000
    R_list = np.arange(0.33, 0.4, 0.01)
    N = len(R_list)
    approx_list = np.zeros(N, dtype=np.float64)
    lit_list = np.zeros(N, dtype=np.float64)

    for rep in range(N_repeat):
        # print(rep, N_repeat)
        for R_i in range(len(R_list)):
            R = R_list[R_i]
            L = 8

            X = np.zeros((int((L/2)**2), 2), dtype=np.float64)
            n_particles = 0

            tries = 0
            while n_particles < (L/2)**2:
                x = np.random.uniform(-L/4, L/4)%L
                y = np.random.uniform(L/4, 3*L/4)
                for n in range(n_particles):
                    rx = np.abs(x - X[n,0])
                    rx = np.minimum(rx, L-rx)
                    ry = np.abs(y - X[n,1])
                    ry = np.minimum(ry, L-ry)

                    r = np.sqrt(rx*rx + ry*ry)
                    if r < 2*R:
                        tries += 1
                        if tries > 5000:
                            n_particles = 0
                            X = np.zeros((L**2, 2), dtype=np.float64)
                            tries = 0
                        break
                else:
                    X[n_particles,0] = x
                    X[n_particles,1] = y
                    n_particles += 1

            A = np.zeros((n_particles, n_particles), dtype=np.float64)
            for p1 in range(n_particles):
                for p2 in range(n_particles):
                    set_value(A, p1, p2, X, L)

            A_diag = np.diag(np.diagonal(A))
            delta_u = np.full(n_particles, 1)

            F_exact = np.linalg.inv(A)@delta_u
            F_approx = np.linalg.inv(A_diag)@delta_u
            F_lit = delta_u

            approx_list[R_i] += np.sqrt(1/n_particles*np.dot(F_exact-F_approx, F_exact-F_approx))
            lit_list[R_i] += np.sqrt(1/n_particles*np.dot(F_exact-F_lit, F_exact-F_lit))

            print(approx_list[R_i], lit_list[R_i])
        
    plt.plot(F_exact)
    plt.show()

    plt.imshow(A.T)
    plt.show()

    plt.scatter(X[:,0], X[:,1])
    plt.show()

    approx_list/=N_repeat
    lit_list/=N_repeat

    plt.plot(R_list, approx_list, label="approx")
    plt.plot(R_list, lit_list, label="literature")
    plt.scatter(R_list, approx_list)
    plt.scatter(R_list, lit_list)
    plt.xlabel("R")
    plt.ylabel("Error")
    plt.legend()
    plt.ylim(bottom=0)
    plt.show()