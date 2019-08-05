import scipy.spatial as sp
import numpy as np 
import matplotlib.pyplot as plt
import numpy.matlib 
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import voronoi3D as voro
import copy
T_max=30
Nx=15
Ny=15
Nz=15
max_x=max(voro.x)
max_y=max(voro.y)
max_z=voro.Zmax
min_x=min(voro.x)
min_y=min(voro.y)
min_z=voro.Zmin
size_x=max_x-min_x
size_y=max_y-min_y
size_z=max_z-min_z
relaxin=True
coeff_diff=1
if relaxin==True:
    coeff_diff=coeff_diff*2

hx=size_x/Nx
hy=size_y/Ny
hz=size_z/Nz
def diffusion(Time_min,Time_max,res_time_min=[1]):
    dt=(1/(8*coeff_diff))*min(hx,hy,hz)**2
    x=np.linspace(min_x,max_x,Nx)
    y=np.linspace(min_y,max_y,Ny)
    z=np.linspace(min_z,max_z,Nz)
    Ux=dt*coeff_diff/hx**2
    Uy=dt*coeff_diff/hy**2
    Uz=dt*coeff_diff/hz**2
    if len(res_time_min)>1:
        Try=copy.deepcopy(res_time_min)
    if Time_min==0:
        Try=np.zeros((Nx,Ny,Nz))
        Try=Try.tolist()
        for j in range(Ny):
            for k in range(Nz):
                Try[0][j][k]=100
    time=Time_min
    while time<Time_max:
        for i in range(1,Nx-1):
            for j in range(1,Ny-1):
                for k in range(1,Nz-1):
                    Try[i][j][k]=Try[i][j][k]+Ux*(Try[i+1][j][k]-2*Try[i][j][k]+Try[i-1][j][k])+Uy*(Try[i][j+1][k]-2*Try[i][j][k]+Try[i][j-1][k])+Uz*(Try[i][j][k+1]-2*Try[i][j][k]+Try[i][j][k-1])
        for u in range(Nx):
            for p in range(Ny):
                Try[0][u][p]=100
        print('time',time)
        time+=dt
    grid=[]
    for i in x :
        for j in y:
            for k in z:
                grid.append([i,j,k])
    printmode=True
    if printmode==True:
        if Time_max==1 or Time_max==3 or Time_max==5:
            X,Y,Z=np.meshgrid(x,y,z)
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111, projection='3d')
            p=[]
            minn, maxx = 0,2
            for i in range(0,Nx):
                for j in range(0,Ny):
                    for k in range(0,Nz):
                        p.append(Try[i][j][k])
            ax.scatter(X,Y,Z,c=p,cmap="inferno",vmin=minn,vmax=maxx)
            colmap =matplotlib.cm.ScalarMappable(cmap="inferno")
            colmap.set_array(p)
            cbar=fig.colorbar(colmap)
            plt.show()
    return(Try,grid)

    

