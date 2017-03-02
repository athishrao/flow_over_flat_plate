

from math import *
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

alpha = np.deg2rad(-8)
Vinf = 45
w=1.7
sina = sin(alpha)

############################################################################################################

x = [.5*(1-(1.1**(17-i))) for i in range(17)]
for i in range(17,37):
    x.append(.05*(i-17))
for i in range(37,59):
    x.append(1-.5*(1-1.1**(i-37)))
##############################################################################################
j = 0
y = []
while (j < 27):
    y.append(-0.15 + 0.5 * (1 - 1.1 ** (26 - j)))
    j+=1

while (j <= 29):
    y.append(0.015 * (j - 29))
    j+=1
y.append(0)  # j=30
j = 31
while j < 34:
    y.append(0.015 * (j - 30))
    j+=1

while j < 60:
    y.append(0.15 - 0.5 * (1 - 1.1 ** (j - 33)))
    j+=1

############################################################################################################
h = []
h.append(None)
for k in range(1, len(x)-1):
    h.append((x[k + 1] - x[k - 1]) / 2.0)

############################################################################################################
k = []
k.append(None)
for m in range(1, 29):
    k.append((y[m + 1] - y[m - 1]) / 2.0)
k.append(None)  # K=29
k.append(None)  # K=30
for m in range(31, 59):
    k.append((y[m + 1] - y[m - 1]) / 2.0)

k[29] = ((y[31] - y[28]) / 2.0)
k.append(None)

############################################################################################################
# Potential at t-1
grid1 = [[0.0 for a in range(60)] for b in range(59)]  # Potential at t=t-1
grid2 = [[0.0 for a in range(60)] for b in range(59)]  # Potential at t=t
V = [[0 for a in range(2)] for b in range(21)]  # Velocity next to the plates
cp = [[1 for a in range(2)] for b in range(21)]  # Values of cp
q = [0 for a in range(21)]


############################################################################################################
def gridval(grid1, i, j):
    # Leakage Factor
    l = (x[58] - x[i]) / (x[58] - x[37])
    # Going through the vertical interior boundaries first
    # LEFT INTERIOR BOUNDARY
    if i == 1:
        # LEFT BOTTOM CORNER
        if j == 1:
            return ((k[j] ** 2.0) * grid1[i + 1][j] + (h[i] ** 2.0) * grid1[i][j + 1]) / ((h[i] ** 2.0) + (k[j] ** 2.0))
        # POINT OF INTERSECTION OF INTERIOR BOUNDRY AND LOWER PLATE
        elif j == 29:
            return ((k[j] ** 2.0) * grid1[i + 1][j] + (h[i] ** 2.0) * (grid1[i][j + 2] + grid2[i][j - 1])) / (2 * (h[i] ** 2.0) + (k[j] ** 2.0))
        # POINT OF INTERSECTION OF INTERIOR BOUNDRY AND UPPER PLATE
        elif j == 30:
            return grid2[i][j - 1]
        # LEFT TOP CORNER
        elif j == 58:
            return ((k[j] ** 2.0) * grid1[i + 1][j] + (h[i] ** 2.0) * grid1[i][j - 1]) / ((h[i] ** 2.0) + (k[j] ** 2.0))
        # NEXT TO LEFT FAR FIELD BOUNDARY
        else:
            return ((k[j] ** 2.0) * grid1[i + 1][j] + (h[i] ** 2.0) * (grid1[i][j + 1] + grid2[i][j - 1])) / (2 * (h[i] ** 2.0) + (k[j] ** 2.0))
    # RIGHT INTERIOR BOUNDARY
    elif i == 57:
        # RIGHT TOP CORNER
        if j == 1:
            return ((k[j] ** 2.0) * grid1[i - 1][j] + (h[i] ** 2.0) * grid1[i][j + 1]) / ((h[i] ** 2.0) + (k[j] ** 2.0))
        # MODIFIED POINT ON INTERIOR BOUNDARY
        elif j == 28:
            return ((k[j] ** 2.0) * (grid2[i - 1][j]) + (h[i] ** 2.0) * (grid1[i][j - 1] - k[j] * l * Vinf * sina)) / ((h[i] ** 2.0) + (k[j] ** 2.0))
        # POINT OF INTERSECTION OF INTERIOR BOUNDRY AND LOWER PLATE
        elif j == 29:
            return grid2[i][j - 1] - k[j - 1] * l * Vinf * sina
        # POINT OF INTERSECTION OF INTERIOR BOUNDRY AND UPPER PLATE
        elif j == 30:
            return grid2[i][j + 1] + k[j + 1] * l * Vinf * sina
        # MODIFIED POINT ON INTERIOR BOUNDARY
        elif j == 31:
            return ((k[j] ** 2.0) * (grid2[i - 1][j]) + (h[i] ** 2.0) * (grid1[i][j + 1] + k[j] * l * Vinf * sina)) / ((h[i] ** 2.0) + (k[j] ** 2.0))
        # RIGHT BELOW CORNER
        elif j == 58:
            return ((k[j] ** 2.0) * grid1[i - 1][j] + (h[i] ** 2.0) * grid1[i][j - 1]) / ((h[i] ** 2.0) + (k[j] ** 2.0))
        # NEXT TO RIGHT FAR FIELD BOUNDARY
        else:
            return ((k[j] ** 2.0) * grid2[i - 1][j] + (h[i] ** 2.0) * (grid1[i][j + 1] + grid2[i][j - 1])) / (2 * (h[i] ** 2.0) + (k[j] ** 2.0))
    # ALL OTHER SPECIAL CONDITIONS
    else:
        # Top Interior Boundary
        if j == 1:
            return ((k[j] ** 2.0) * (grid1[i + 1][j] + grid2[i - 1][j]) + (h[i] ** 2.0) * grid1[i][j + 1]) / ((h[i] ** 2.0) + 2 * (k[j] ** 2.0))
        # LAYER JUST BELOW LOWER LAYER OF PLATE
        elif j == 28:
            # DOWNSTEAM FROM PLATE
            if i < 17:
                return ((k[j] ** 2.0) * (grid1[i + 1][j] + grid2[i - 1][j]) + (h[i] ** 2.0) * (grid1[i][j + 1] + grid2[i][j - 1])) / (2 * (h[i] ** 2 + k[j] ** 2.0))
            # ON THE PLATE
            elif i <= 37:
                return ((k[j] ** 2.0) * (grid1[i + 1][j] + grid2[i - 1][j]) + (h[i] ** 2.0) * (grid1[i][j - 1] - k[j] * Vinf * sina)) / ((h[i] ** 2.0) + 2 * (k[j] ** 2.0))
            # UPSTREAM TO PLATE
            else:
                return ((k[j] ** 2.0) * (grid1[i + 1][j] + grid2[i - 1][j]) + (h[i] ** 2.0) * (grid1[i][j - 1] - k[j] * l * Vinf * sina)) / ((h[i] ** 2.0) + 2 * (k[j] ** 2.0))
        # UPPER LAYER OF THE PLATE
        elif j == 29:
            # DOWNSTEAM FROM PLATE
            if i < 17:
                return ((k[j] ** 2.0) * (grid1[i + 1][j] + grid2[i - 1][j]) + (h[i] ** 2.0) * (grid1[i][j + 2] + grid2[i][j - 1])) / (2 * (h[i] ** 2 + k[j] ** 2.0))
            # ON THE PLATE
            elif i <= 37:
                return grid2[i][j - 1] - k[j - 1] * Vinf * sina
            # UPSTREAM TO PLATE
            else:
                return grid2[i][j - 1] - k[j - 1] * l * Vinf * sina
        # LOWER LAYER OF THE PLATE
        elif j == 30:
            # DOWNSTEAM FROM PLATE
            if i < 17:
                return grid2[i][j - 1]
            # ON THE PLATE
            elif i <= 37:
                return grid2[i][j + 1] + k[j + 1] * Vinf * sina
            # UPSTREAM TO PLATE
            else:
                return grid2[i][j + 1] + k[j + 1] * l * Vinf * sina
        # LAYER JUST ABOVE UPPER LAYER OF PLATE
        elif j == 31:
            # DOWNSTEAM FROM PLATE
            if i < 17:
                return ((k[j] ** 2.0) * (grid1[i + 1][j] + grid2[i - 1][j]) + (h[i] ** 2.0) * (grid1[i][j + 1] + grid2[i][j - 1])) / (2 * (h[i] ** 2 + k[j] ** 2.0))
            # ON THE PLATE
            elif i <= 37:


                return ((k[j] ** 2.0) * (grid1[i + 1][j] + grid2[i - 1][j]) + (h[i] ** 2.0) * (grid1[i][j + 1] + k[j] * Vinf * sina)) / ((h[i] ** 2.0) + 2 * (k[j] ** 2.0))
            # UPSTREAM TO PLATE
            else:
                return ((k[j] ** 2.0) * (grid1[i + 1][j] + grid2[i - 1][j]) + (h[i] ** 2.0) * (grid1[i][j + 1] + k[j] * l * Vinf * sina)) / ((h[i] ** 2.0) + 2 * (k[j] ** 2.0))
        # Below Interior Boundary
        elif j == 58:
            return ((k[j] ** 2.0) * (grid1[i + 1][j] + grid2[i - 1][j]) + (h[i] ** 2.0) * grid2[i][j - 1]) / ((h[i] ** 2.0) + 2 * (k[j] ** 2.0))
        # ALL THE INTERIOR POINTS
        else:
            return ((k[j] ** 2.0) * (grid1[i + 1][j] + grid2[i - 1][j]) + (h[i] ** 2.0) * (grid1[i][j + 1] + grid2[i][j - 1])) / (2 * (h[i] ** 2 + k[j] ** 2.0))


def intrapolated(grid2, i, j, w):
    grid2[i][j] = (1 - w) * grid1[i][j] + w * gridval(grid2, i, j)


def conv():
    e = 0
    for i in range(1, 58):
        for j in range(1, 59):
            if abs(grid1[i][j] - grid2[i][j]) > e:
                e = abs(grid1[i][j] - grid2[i][j])
    return e


############################################################################################################
t = 0
b = True
e = 1
eps = 0.001
print 'sasaa',x[17],x[37]
while b:
    for i in range(1, 58):
        for j in range(1, 59):
            grid1[i][j] = (grid2[i][j])
    for i in range(1, 58):
        for j in range(1, 59):
            intrapolated(grid2, i, j, w)
    for i in range(17, 58):
        intrapolated(grid2, i, 30, w)
    grid2[0][0] = grid2[1][1]
    grid2[0][59] = grid2[1][58]
    grid2[58][0] = grid2[57][1]
    grid2[58][59] = grid2[57][58]
    for j in range(1, 59):
        grid2[0][j] = grid2[1][j]
        grid2[58][j] = grid2[57][j]
    for i in range(1, 58):
        grid2[i][0] = grid2[i][1]
        grid2[i][59] = grid2[i][58]
    t += 1
    e = conv()
    if e > eps:
        b = True
    else:
        b = False
    if t>700:
        break
    print ("t= ", t)

print 'sasaa', x[17], x[37]


############################################################################################################
def velocity_cp(a,grid):
    v = [[0.0 for i in range(21)] for j in range(2)]
    cp =[[0.0 for i in range(21)] for j in range(2)]
    v[0][0] = ((grid[18][30] - grid[16][30]) / .1) + (Vinf * cos(alpha))
    v[1][0] = ((grid[18][29] - grid[16][29]) / .1) + (Vinf * cos(alpha))
    v[0][20] = ((grid[38][30] - grid[36][30]) / .1) + (Vinf * cos(alpha))
    v[1][20] = ((grid[38][29] - grid[36][29]) / .1) + (Vinf * cos(alpha))

    cp[0][0] = 1-((v[0][0]/45.0)**2)
    cp[1][0]= 1-((v[1][0]/45.0)**2)
    cp[0][20]= 1-((v[0][20]/45.0)**2)
    cp[1][20]= 1-((v[1][20]/45.0)**2)
    for i in range(1, 20):
        v[1][i] = ((grid[i + 18][29] - grid[i+16][29]) / (.1)) + (Vinf * cos(alpha))
        cp[1][i] = (1-((v[1][i]/45.0)**2))
        v[0][i] = ((grid[i+18][30] - grid[i+16][30]) / (.1)) + (Vinf * cos(alpha))
        cp[0][i] = (1 - ((v[0][i] / 45.0) ** 2))
    if a==0:
        return v
    else:
        return cp


def delta_cp():
    cp = velocity_cp(1,grid2)
    k=[]
    for i in range(21):
        k.append(cp[1][i]-cp[0][i])
    return k

def coefflift():
    diff = delta_cp()
    u = (1/21.0) * sum(diff)
    return u
error = abs(coefflift()-(2*pi*alpha))
print 'error precentage',abs(error)/abs(2*pi*alpha)
print('FINAL!!!!!!!',error)
#gridblu=[list(x) for x in zip(*grid2)]
############################################################################################################
v =velocity_cp(0,grid2)
cp = velocity_cp(1,grid2)
plt.plot(x[17:38],v[0],'-r')
plt.plot(x[17:38],v[1],'-b')
plt.ylabel('Upper and Lower Velocities')
plt.xlabel('x[i]')
plt.show()
plt.plot(x[17:38],cp[0],'-r')
plt.plot(x[17:38],cp[1],'-b')
plt.ylabel('Upper and lower Pressures')
plt.xlabel('x[i]')
plt.show()
plt.plot(x[17:38],delta_cp(),'-g')
plt.ylabel('Delta Cp')
plt.xlabel('x[i]')
plt.show()
gridblu=[list(x) for x in zip(*grid2)]
fig = plt.figure()
axes=Axes3D(fig)
X,Y = np.meshgrid(x,y)
axes.plot_surface(X,Y,gridblu,rstride=1,cstride=1,cmap=plt.cm.spectral,linewidth=0,antialiased=False)
plt.xlabel('x[i]')
plt.ylabel('y[i]')
axes.set_zlabel('Phi')
plt.show()
############################################################################################################
