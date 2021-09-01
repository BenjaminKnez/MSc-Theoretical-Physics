import random

import numpy as np
import math


L0 = 1.0
dL = 0.01
dt = 0.01
D0 = 0.01
T0 = 100000
iter = 300

c1 = 0.1
c2 = 0.5

ATP0 = 100

length = np.zeros(T0)   #large size

for i in range(1, iter):

    l = np.random.exponential(L0, None)

    ATP_count = 0

    for t in range(0, T0):

        b1 = -l / 2
        b2 = l / 2

        x = random.uniform(b1, b2)

        for j in range(0, int((x-b1)/dL)):

            if random.uniform(0, 1)<2*c1*dL*dt:
                b1 = x - j*dL
                break

        for k in range(0, int((b2-x)/dL)):

            if random.uniform(0, 1)<2*c1*dL*dt:
                b2 = x + k*dL
                break

        L1 = np.random.exponential(l, None)
        if random.uniform(0, 1)<(c2*dt*(ATP0-ATP_count)/ATP0) and ATP_count < ATP0:
            b1 = b1 - L1
            ATP_count = ATP_count + 1

        L2 = np.random.exponential(l, None)
        if random.uniform(0, 1)<(c2*dt*(ATP0-ATP_count)/ATP0) and ATP_count < ATP0:
            b2 = b2 + L2
            ATP_count = ATP_count + 1

        l = b2 - b1

        if l <= 0:
            break

        length[t] = length[t] + l/iter

print(length)

np.savetxt("Living_ATP.txt", length)


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


# xvalues = np.array([0.5,1,2,5,10,20])
# yvalues = np.zeros(6)
# for i in range(6):
#     yvalues[i] = Rep_living(xvalues[i], 1.0)
#
# plt.plot(xvalues, yvalues, 'bo', linestyle="None")
# plt.yscale('log')
# plt.xscale('log')
# plt.title("Dependence of viscosity on breakage rate (xi)")
# plt.xlabel("Adimensional breakage rate")
# plt.ylabel("Adimensional viscosity")
# plt.ylim(1, 1000)
# plt.show()
