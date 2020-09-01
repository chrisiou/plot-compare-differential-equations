# Finite differences for a Dirichlet bvp

# Suppose y(x)=cos(x) the solution of the boundary-value problem: -y''(x)+y(x)=2cos(x) forall x in [0,π] and
# with y(0)=1, y(π)=-1. Solve the bvp using Ν=50 and plot the approximate and the actual solutions.
# Furthermore, for N = [100,200,400] find the approach errors.

import numpy as np
import matplotlib.pyplot as plt

def y(x):
    return np.cos(x)

def q(x):
    return 1

def f(x):
    return 2*np.cos(x)

def fdm(N, yA, yB):
    h = np.pi/N
    t = np.linspace(0, np.pi, N+1)

    # definition of array F, (AU=F)
    F = f(t[1:-1])*(h**2)
    F[0] += yA
    F[-1] += yB

    Q = np.ones(N-1)*q(t[1:-2])*(2+h**2)
    Z = np.ones(N-2)*(-1)

    # AU = F => (Z+h(^2)Q)U = F(h^2)
    A = np.diag(Q,0) + np.diag(Z, 1) + np.diag(Z, -1)

    U = np.linalg.solve(A,F)
    U = np.insert(U, 0, yA)
    U = np.append(U, yB)
    
    return U
    
def makeGraph(U, N):
    t = np.linspace(0, np.pi, N+1)

    plt.plot(U, t, linewidth = 5, label = "approx sol")
    plt.plot(y(t), t, marker='o', markersize=4, label = "actual sol")

    plt.xlabel('t')
    plt.ylabel('y(t)')
    plt.title("for N = " + str(N))
    plt.legend()
    plt.show()

yA = 1
yB = -1
N = 50

U = fdm(N, yA, yB)

# Errors
N2 = [100, 200, 400]
for n in N2:
    t = np.linspace(0, np.pi, n+1)
    print("for N =", n, " error is: ", np.linalg.norm(fdm(n, yA, yB) - y(t), np.inf))

# convergence analysis: (np.log(E[i+1]/E[i]) / np.log( N[i] / N[i+1])), for E a list of errors for different N

makeGraph(U, N)
