import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d

def isConvergent(m1, m2, e):
    for i in range(len(m1)):
        for j in range(len(m1[0])):
            if abs(m1[i][j]-m2[i][j]) > e:
                return False
    return True
'''Used for comparrison'''
itJ = []
itGS = []
itGSOR = []
cCriteria = []
'''Used for comparrison'''
for o in range(1, 16): #Solving for different convergence criterias
    
    xGrid = np.linspace(0, np.pi, 11)
    yGrid = np.linspace(0, np.pi/2, 21)
    h = xGrid[1] - xGrid[0]
    k = yGrid[1] - yGrid[0]
    beta = h/k
    
    u0 = np.zeros((len(xGrid), len(yGrid)))

    cCriteria.append(1*(10**(-o)))
    
    '''{Applying boundary conditions'''
    u0[0] = [np.cos(i) for i in yGrid]
    u0[-1] = [-np.cos(i) for i in yGrid]
    for i in range(len(u0)):
        u0[i][0] = np.cos(xGrid[i])
        u0[i][-1] = 0
    '''Applying boundary conditions}'''
    
    '''{Jacobi'''
    u = np.copy(u0)
    jac = np.copy(u0)
    it = 0
    while True:
        it+=1
        u = np.copy(jac)  #Values from old iteration
        for i in range(1, len(u)-1):
            for j in range(1, len(u[0])-1):
                jac[i][j] = (u[i+1][j] + u[i-1][j] + (beta**2)*(u[i][j+1] + u[i][j-1] + 
                 (h**2)*(np.cos(xGrid[i] + yGrid[j]) + np.cos(xGrid[i] - yGrid[j]))))/(2+2*(beta**2)) 
        if isConvergent(u, jac, 1*(10**(-o))):
            break
    itJ.append(it)
    '''Jacobi}'''
    
    '''{Gauss-Seidel'''
    gs = np.copy(u0)
    it = 0
    while True:
        it+=1
        gsoi = np.copy(gs)  #Values from old iteration
        for i in range(1, len(u)-1):
            for j in range(1, len(u[0])-1):
                gs[i][j] = (gs[i+1][j] + gs[i-1][j] + (beta**2)*(gs[i][j+1] + gs[i][j-1] + 
                 (h**2)*(np.cos(xGrid[i] + yGrid[j]) + np.cos(xGrid[i] - yGrid[j]))))/(2+2*(beta**2)) 
        if isConvergent(gsoi, gs, 1*(10**(-o))):
            break
    itGS.append(it)
    '''Gauss-Seidel}'''
        
    '''{GSSOR'''
    gsor = np.copy(u0)
    it = 0
    while True:
        it+=1
        gsoroi = np.copy(gsor) #Values from old iteration
        for i in range(1, len(u)-1):
            for j in range(1, len(u[0])-1):
                gsor[i][j] = -0.25*gsoroi[i][j] + 1.25*(gsor[i+1][j] + gsor[i-1][j] + (beta**2)*(gsor[i][j+1] + 
                 gsor[i][j-1] + (h**2)*(np.cos(xGrid[i] + yGrid[j]) + np.cos(xGrid[i] - yGrid[j]))))/(2+2*(beta**2)) 
        if isConvergent(gsoroi, gsor, 1*(10**(-o))):
            break
    itGSOR.append(it)
    '''GSSOR}'''
        
'''{Plot'''
fig2 = plt.figure()
plt.plot(cCriteria, itJ, label='Jacobi')
plt.plot(cCriteria, itGS, label='Gauss-Seidel')
plt.plot(cCriteria, itGSOR, label = 'GSSOR')
plt.legend(loc='best')
plt.xlabel('Convergence criteria')
plt.ylabel('Iteration count')
plt.xscale('log')
plt.show()
fig2.set_size_inches(16,8)
fig2.savefig('hw3q1b', dpi=100)

X, Y = np.meshgrid(xGrid, yGrid)
fig = plt.figure()
ax = fig.gca(projection='3d')
plt.xlabel('x')
plt.ylabel('y')
fig.set_size_inches(16,8)
surf = ax.plot_surface(X, Y ,np.transpose(u), cmap=cm.coolwarm)
fig.savefig('hw3q1.png', dpi=100)
'''Plot}'''