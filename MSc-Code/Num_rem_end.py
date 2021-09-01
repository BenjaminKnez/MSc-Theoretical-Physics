#start with average length, and at each time appropriately scale the diffusion coefficeint (length) of ALL ends

import sys

import random

import numpy as np
import math
import matplotlib.pyplot as plt

#This code works only when initial chains are longer than the dimensions of the box, for the approximation of random ends positions to be sensible

#usual constants

L0 = 150.0
dL = 0.01
dt = 0.01
D0 = 0.1
c1 = 0.1


N0 = 40    #number of initial chains
A = 5.0     #size of the box
T0 = 500000 #time of observation
X = 0.3    #distance between ends for ligation to occur, size of ligation enzyme?


#matrices are t - k - i

ends = np.zeros((T0, 3, 2*N0))

surv = np.zeros(T0-1)


L = L0

n = 0

for j in range(2*N0):

    for k in range(3):

        ends[0][k][j] = random.uniform(0, A)


for t in range(0, T0-1):


    surv[t] = 2*N0 - n

    ends[t] = ends[t][np.arange(ends[t].shape[0])[:, None], (ends[t] == 0).argsort(1, kind="mergesort")]

    L = L0*(N0*1.0/(surv[t]/2.0))

    for i in range(0, int(surv[t])):

        dir = random.choice([0, 1, 2])

        ends[t+1][dir][i] = ends[t][dir][i] + random.choice([-1, 1]) * math.sqrt(2.0*D0*dt/L)
        if ends[t+1][dir][i] > A:
            ends[t+1][dir][i] = ends[t][dir][i] - 2 * math.sqrt(2.0 * D0 * dt / L)
        elif ends[t + 1][dir][i] < 0:
            ends[t + 1][dir][i] = ends[t][dir][i] + 2 * math.sqrt(2.0 * D0 * dt / L)

        ends[t + 1][dir-1][i] = ends[t][dir-1][i]
        ends[t + 1][dir-2][i] = ends[t][dir-2][i]

    for m in range(int(surv[t])):

        if ends[t+1][0][m] > 0 and ends[t+1][1][m] > 0 and ends[t+1][2][m] > 0:

            for j in range(m+1, int(surv[t])):

                if ends[t + 1][0][j] > 0 and ends[t + 1][1][j] > 0 and ends[t + 1][2][j] > 0 and math.sqrt((ends[t + 1][0][j] - ends[t + 1][0][m]) ** 2 + (ends[t + 1][1][j] - ends[t + 1][1][m]) ** 2 + (ends[t + 1][2][j] - ends[t + 1][2][m]) ** 2) < X and random.uniform(0, 1) < c1 * dL * dt:

                    n = n + 2

                    for k in range(3):
                        ends[t+1][k][m] = 0
                        ends[t + 1][k][j] = 0   #random number to put them out of the run

                    break



print(surv)

np.savetxt("Num_" + str(sys.argv[1]) +".txt", surv)


# xaxis = np.arange(0, T0-1)

# plt.plot(xaxis, surv)
# plt.xlabel("Time [t]")
# plt.ylabel("Number of remaining ends")
# plt.ylim(0, 2*N0)
# plt.show()
