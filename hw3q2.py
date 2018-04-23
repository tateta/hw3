import numpy as np
import matplotlib.pyplot as plt
def converges(m1, m2, e):
    for i in range(len(m1)):
        for j in range(len(m1[0])):
            if abs((m1-m2)[i][j]) > e:
                return False
    return True

h = 0.025
re = 500

def withoutUpwindDownwind(h, re):
    xGrid = np.linspace(0, 1, int(1/h+1))
    yGrid = np.linspace(0, 2, int(2/h+1))
    u0 = np.zeros((len(xGrid), len(yGrid)))
    for i in range(len(xGrid)):
        u0[i][-1] = 1
    u = np.copy(u0)
    a = 2/(re*h)
    it = 0
    upi = np.copy(u)
    while True:
        it+=1
        upi = np.copy(u)
        for i in range(1, len(u)-1):
            for j in range(1, len(u[0])-1):
                u[i][j] = -u[i][j]*(u[i+1][j] - u[i-1][j])/(4*a) + 0.25*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1])
        if it%50==0:
            print(it, np.sum(np.abs(upi-u)))
        if converges(u, upi, 1e-8):
            break
    return u

def withUpwindDownwind(h, re):
    xGrid = np.linspace(0, 1, int(1/h+1))
    yGrid = np.linspace(0, 2, int(2/h+1))
    u0 = np.zeros((len(xGrid), len(yGrid)))
    for i in range(len(xGrid)):
        u0[i][-1] = 1
    u = np.copy(u0)
    a = 2/(re*h)
    it = 0
    upi = np.copy(u)
    while True:
        it+=1
        upi = np.copy(u)
        for i in range(1, len(u)-1):
            for j in range(1, len(u[0])-1):
                if True:
                    u[i][j] = -u[i][j]*(u[i][j] - u[i-1][j])/(4*a) + 0.25*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1])
                else:
                    u[i][j] = -u[i][j]*(u[i+1][j] - u[i][j])/(4*a) + 0.25*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1])
        if converges(u, upi, 1e-8):
            break
        if it%50==0:
            print('ud', it, np.sum(np.abs(upi-u)))
    return u

u = withoutUpwindDownwind(h, re)
fig = plt.figure()
CS=plt.contour(np.linspace(0, 1, int(1/h+1)),np.linspace(0, 2, int(2/h+1)),np.transpose(u))
plt.clabel(CS, inline=1)
fig.set_size_inches(6,12)
fig.savefig('hw3q2R5.png', dpi=100)
plt.show()