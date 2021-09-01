import sys

import random

import numpy as np
import math
import matplotlib.pyplot as plt



def Rep_living(ksi, L0):

    dL = 0.01
    dt = 0.01
    D0 = 0.01
    T0 = 30000
    iter = 500


    c1 = D0/(L0**4 * ksi)

    global surv
    surv = np.zeros(T0)   #large size

    global length
    length = np.zeros(T0)

    for i in range(1, iter):

        L = np.random.exponential(L0, None)
        b1 = -L/2
        b2 = L/2
        x = random.uniform(b1, b2)

        for t in range(0, T0):

            x = x + random.choice([-1, 1]) * math.sqrt(2*D0*dt/L)

            if x<b1:
                break

            if x>b2:
                break

            for j in range(0, int((x-b1)/dL)):

                if random.uniform(0, 1)<2*c1*dL*dt:
                    b1 = x - j*dL
                    break

            for k in range(0, int((b2-x)/dL)):

                if random.uniform(0, 1)<2*c1*dL*dt:
                    b2 = x + k*dL
                    break

            L1 = np.random.exponential(L0, None)
            if random.uniform(0, 1)<c1*dt:
                b1 = b1 - L1

            L2 = np.random.exponential(L0, None)
            if random.uniform(0, 1)<c1*dt:
                b2 = b2 + L2

            L = b2 - b1

            surv[t] = surv[t] + 1.0/iter

            length[t] = length[t] + L / iter

    print(surv)

    print(length)

    global visc

    visc = 0

    for l in range(0, T0):

        visc = visc + (surv[l]) * dt

    print(visc)

    return visc

Rep_living(float(sys.argv[1]), 1.0)

np.savetxt("Surv_Ksi_" + str(sys.argv[1]) + ".txt", surv)
np.savetxt("Visc_Ksi_" + str(sys.argv[1]) + ".txt", visc)

# xaxis = np.arange(0, T0)

# Rep_living(100000, 1.0)  #non breaking
# plt.plot(xaxis, surv, 'k')
#
# Rep_living(100.0, 1.0)
# plt.plot(xaxis, surv, 'b')
#
# Rep_living(10.0, 1.0)
# plt.plot(xaxis, surv, 'g')
#
# Rep_living(1.0, 1.0)
# plt.plot(xaxis, surv, 'r')
#
# Rep_living(0.1, 1.0)
# plt.plot(xaxis, surv, 'c')
#
# Rep_living(0.01, 1.0)
# plt.plot(xaxis, surv, 'm')




# xvalues = np.array([0.001,0.005,0.01,0.05,0.1,0.5,1.0,5.0,10.0,50.0,100.0])
# yvalues = np.zeros(11)
#
# for i in range(11):
#     yvalues[i] = Rep_living(xvalues[i], 1.0)
#
#
#
# plt.plot(xvalues, yvalues, 'bo', linestyle="None")
# plt.yscale('log')
# plt.xscale('log')
# plt.title("Dependence of viscosity on breakage rate (xi)")
# plt.xlabel("Adimensional breakage rate")
# plt.ylabel("Adimensional viscosity")
# plt.ylim(1, 1000)
# plt.show()
