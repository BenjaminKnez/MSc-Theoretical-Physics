import sys
import random

import numpy as np
import math
from scipy.optimize import fsolve

def Rep_fusion(phi, Ta):

    n0 = 1000.0
    ATP0 = 500.0

    L0 = 1.0
    dL = 0.01
    dt = 0.01
    D0 = 0.01
    iter = 2000
    T0 = 1000000

    ksi = 1.0/phi
    c1 = D0/(L0**4 * ksi)

    Td = L0**3 / D0

    A = (phi/Td)*((ATP0-n0)/ATP0)
    B = (phi/(Td*ATP0))

    def f(n):
        return (B / A ** 2) * math.log(abs(n / (B * n + A))) + (1 / A * n) - (B / A ** 2) * math.log(abs(n0 / (B * n0 + A))) - (  1 / A * n0) - Ta


    a=fsolve(f,100)


    global surv
    surv = np.zeros(T0)   #large size

    for i in range(1, iter):

        L = np.random.exponential((L0*n0)/a[0], None)   #discrete n, should depend on concentration of enzymes and ends
        b1 = -L/2
        b2 = L/2
        x = random.uniform(b1, b2)

        for t in range(0, T0):

            x = x + random.choice([-1, 1]) * math.sqrt(2*D0*dt/L)

            if x<b1:
                break

            if x>b2:
                break

            # L1 = np.random.exponential(L0*math.exp(n), None)
            # if random.uniform(0, 1) < c1 * dL * dt:
            #     b1 = b1 - L1
            #
            # L2 = np.random.exponential(L0*math.exp(n), None)
            # if random.uniform(0, 1) < c1 * dL * dt:
            #     b2 = b2 + L2
            #
            #
            # L = b2 - b1

            surv[t] = surv[t] + 1.0/iter

    visc = 0

    for l in range(0, T0):

        visc = visc + (surv[l]) * dt

    print(visc)

    return visc



Rep_fusion(float(sys.argv[1]), float(sys.argv[2]))
