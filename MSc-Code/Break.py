import sys
import random

import numpy as np
import math


def Rep_breaking(chi, Ta):

    L0 = 1.0
    dL = 0.01
    dt = 0.01
    D0 = 0.01
    iter = 500
    T0 = 10000
    #Td = 2s, Tbreak = 12-36min

    ksi = 1.0/chi
    c1 = D0/(L0**4 * ksi)

    Td = L0**3 / D0
    n = (chi*Ta)/Td

    global surv
    surv = np.zeros(T0)  #large size

    for i in range(1, iter):

        L = np.random.exponential(L0/(n+1), None)
        b1 = -L/2
        b2 = L/2
        x = random.uniform(b1, b2)

        for t in range(0, T0):

            x = x + random.choice([-1, 1]) * math.sqrt(2*D0*dt/L)

            if x<b1:
                break

            if x>b2:
                break

            # for j in range(0, int((x-b1)/dL)):
            #
            #     if random.uniform(0, 1)<c1*dL*dt:
            #         b1 = x - j*dL
            #         break
            #
            # for k in range(0, int((b2-x)/dL)):
            #
            #     if random.uniform(0, 1)<c1*dL*dt:
            #         b2 = x + k*dL
            #         break
            #
            #
            #
            # L = b2 - b1

            surv[t] = surv[t] + 1.0/iter

    visc = 0

    for l in range(0, T0):

        visc = visc + (surv[l]) * dt

    print(visc)

    return visc



Rep_breaking(float(sys.argv[1]), float(sys.argv[2]))
