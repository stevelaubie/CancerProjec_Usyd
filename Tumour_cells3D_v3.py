from shapely.geometry import LineString, MultiPolygon, MultiPoint, Point
import scipy.spatial as sp
import numpy as np 
import matplotlib.pyplot as plt
import numpy.matlib 
import voronoi3D as V
import Healthy_cells3D as HC
import Diffusion3D_v2 as D
import pyvoro
from pyvoro import voroplusplus
import pickle
import random
import os
import cv2
from cv2 import VideoWriter, VideoWriter_fourcc, imread, resize
import copy

Allcells={}



class All_cells():
    
    def __init__(self,cellid,pos,celltype=0,collagen_value=0,spring=2.8,birthtime=None,shift=0):
        self.cellid=cellid
        self.pos=pos
        self.celltype=celltype
        self.collagen_value=collagen_value
        self.spring=spring
        #for tumour cells only
        self.birthtime=birthtime
        self.shift=shift
        
    def distance(self,other):
        return np.sqrt((self.pos[0]-other.pos[0])**2+(self.pos[1]-other.pos[1])**2+(self.pos[2]-other.pos[2])**2)
    
    def changer_position(self,i,j,k):
        self.pos[0]=i
        self.pos[1]=j
        self.pos[2]=k
    
    def abscisse(self):
        return self.pos[0]

    def ordonnee(self):
        return self.pos[1]
    
    def profondeur(self):
        return self.pos[2]

    #def type_of_cell(self):
     #   if self.celltype==0:
           # print("This is a healthy cell")
      #  else:
           # print("This is a tumour cell")
        
    def age(self,t):
        if self.celltype==1:
            return t-self.birthtime
    
    def find_cellids_of_neighbors(self,k1):
        #returns a list of all the cell_ids of the neighbors.
        neighbors=[]
        point_index=self.cellid-self.shift
        cell=k1[point_index]
        for faces in cell["faces"]:
            if faces['adjacent_cell']>-1:
                u=0
#                for k in range(len(Allcells)):
#                    if Allcells[k].cellid==faces['adjacent_cell']:
#                        u=Allcells[k].shift
                neighbors.append(faces['adjacent_cell'])
        return list(set(neighbors))
    
    def cos_theta(self,other):
        x=(other.pos[0]-self.pos[0])/self.distance(other)
        return x
    
    def sin_theta(self,other):
        y=(other.pos[1]-self.pos[1])/self.distance(other)
        return y
    def Z_theta(self,other):
        z=(other.pos[2]-self.pos[2])/self.distance(other)
        return z



next_cellid=len(V.coord)

def getnext_cellid():
    global next_cellid
    r = next_cellid
    next_cellid += 1
    return r 

for i in range(len(V.coord)):
    Allcells[i]=All_cells(i,[V.coord[i][0],V.coord[i][1],V.coord[i][2]])

'''for i in range(len(E.coordonnees_bis)):
    print(i,E.coordonnees_bis[i])'''



#It's the function that will influence the tumour cells to go one direction or another.

class virus():
    def __init__(self,cellid,nb,shift):
        self.cellid=cellid
        self.nb=nb
        self.shift=shift


def f(x):return 80*x+1 


#rgrowth is the tumor growth parameter
rgrowth=1.5
s=0

def add_C_cell(c):
        global next_cellid
        keyss=0
        if c.cellid not in Allcells.keys():
            for keys in Allcells.keys():
                if keys>keyss:
                    keyss=keys
#            Allcells[c.cellid-c.shift]=c
            Allcells[keyss+1]=c
        if c.cellid>=next_cellid:
           # print(type(next_cellid))
            #print(type(c.cellid))
            next_cellid=c.cellid+1


def closest_element(l,a):
    s=abs(l[0]-a)
    d,g=0,0
    for i in l:
        if abs(i-a)<s:
            s,g=abs(i-a),i-a
            d+=1
    if g<0 :return d+1
    else: return d
    
def friendly_neighbours(a,p,k):
    friendly_list=[]
    id=[]
    for i in range(1,len(k)):
        
         if k[i].distance(a)<p :
             if k[i].celltype==1:
                 friendly_list.append(k[i])
                 id.append(k[i].cellid)
    return(id)

def distance3D(a,b):
    return(np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2))

def closest_3Dposition(PositionCell,Grid):
    mindistance=1000
    for i in range (len(Grid)):
        if distance3D(PositionCell,Grid[i])<mindistance:
            mindistance=distance3D(PositionCell,Grid[i])
            index=i
    return(index)
        
    
def list_of_intervals(l):
        s=l[0]
        CA_collagen_intervals=[]
        CA_collagen_intervals.append(s)
        for i in range (1,len(l)):
            s+=l[i]
            CA_collagen_intervals.append(s)
        return CA_collagen_intervals



#print(len(E.coordonnees_bis))
epsilon=0.8
dt=1


Allcells[450].celltype=1
Allcells[450].birthtime=-1
first_cancer_cell=Allcells[450]

tri=V.coord
data=[]
data.append(V.k)
last_lysis=91919
counterr=0

def f(x,y,z):
    return np.sqrt(x**2+y**2+z**2)

def distance_deux_points(x,y,z,a,b,c):
    return np.sqrt((x-a)**2+(y-b)**2+(z-c)**2)



def g(x):return 80*x+1


ovarianshape=False
if ovarianshape==False:
#    K=[[20.784609690826528,-9.0],[20.784609690826528,-6.0],[20.784609690826528,-3.0],[20.784609690826528,0.0],[20.784609690826528,3.0],[20.784609690826528,6.0],[20.784609690826528,9.0],[20.784609690826528, 12.0],[18.186533479473212 ,-8.0],
#    [18.186533479473212,-5.0],[18.186533479473212,-2.0],[18.186533479473212,1.0],[18.186533479473212,4.0],[18.186533479473212,7.0],[18.186533479473212 ,10],[20.784609690826528,51],[18.186533479473212,55.0],[18.186533479473212 ,52.0],[12.990381056766578,64.0],[23.38268590217984, 55.0],
#    [23.38268590217984, -8.0],[23.38268590217984, -5.0],[23.38268590217984,-2.0],[23.38268590217984, 1.0],[23.38268590217984 , 4.0],[23.38268590217984 , 7.0],[20.784609690826528, 48.0 ],[23.38268590217984, 52.0],
#    [20.784609690826528,51.0],[23.38268590217984 ,10.0],[20.784609690826528,54.0],[25.980762113533157 ,57.0],[25.980762113533157, 60.0],[25.980762113533157,63.0],[20.784609690826528,57.0],[20.784609690826528,60.0],[20.784609690826528,63.0],[18.186533479473212,58.0],[18.186533479473212, 61.0],[18.186533479473212 ,64.0],
#    [23.38268590217984 ,58],[23.38268590217984, 61.0],[23.38268590217984,64.0],[15.588457268119894, 60.0],[15.588457268119894,63.0],[36.373066958946424,21.0],
#    [36.373066958946424,24.0],[36.373066958946424,27.0],[36.373066958946424,30.0],[36.373066958946424,33.0],[36.373066958946424,36.0],[28.578838324886473, 64.0],
#    [36.373066958946424,39.0],[38.97114317029974 ,22.0],[38.97114317029974, 25.0],[38.97114317029974 ,28.0],[38.97114317029974 ,31.0],[38.97114317029974 ,34.0],[38.97114317029974 ,37.0],[38.97114317029974 ,40.0],[41.569219381653056, 21.0], [44.167295593006365 ,40.0],[49.363448015713,40.0],[41.569219381653056,24.0],[41.569219381653056,27.0],[41.569219381653056,30.0],
#    [41.569219381653056,33.0],[41.569219381653056,36.0],[41.569219381653056,39.0],[44.167295593006365,22.0],[44.167295593006365,
#    .0],[44.167295593006365,28.0],[44.167295593006365,31.0],[44.167295593006365,34.0],[44.167295593006365,37.0],[46.76537180435968,21.0],
#    [46.76537180435968,24.0],[46.76537180435968,27.0],[46.76537180435968,30.0],[46.76537180435968,33.0],[46.76537180435968,36.0],[46.76537180435968,39.0],[49.363448015713,22.0],[49.363448015713,25.0],[49.363448015713,28.0],[49.363448015713,31.0],
#    [49.363448015713,34.0],[49.363448015713,37.0],[49.363448015713 ,22.0],[49.363448015713 ,25.0],[49.363448015713 ,28.0],[49.363448015713 ,31.0],[49.363448015713,34.0],[49.363448015713,37.0],[49.363448015713,40.0],[51.96152422706631,21.0],
#    [51.96152422706631 ,24.0],[51.96152422706631 ,27.0],[51.96152422706631 ,30.0],[51.96152422706631,33.0],[51.96152422706631,36.0],[51.96152422706631 ,39.0], [0.0, 15.0], [0.0, 18.0], [0.0, 21.0],
#    [0.0, 24.0], [0.0, 27.0], [0.0, 30.0], [0.0, 33.0], [2.598076211353316, 16.0], [2.598076211353316, 19.0], [2.598076211353316, 22.0], [2.598076211353316, 25.0],
#    [2.598076211353316, 28.0], [2.598076211353316, 31.0], [5.196152422706632, 18.0], [5.196152422706632, 21.0], [5.196152422706632, 24.0], [5.196152422706632, 27.0],
#    [5.196152422706632, 30.0], [7.794228634059947, 19.0], [7.794228634059947, 22.0], [7.794228634059947, 25.0], [7.794228634059947, 28.0],[10.392304845413264, 21.0],[10.392304845413264, 24.0],[10.392304845413264, 27.0],[12.990381056766578, 22.0],
#    [12.990381056766578,25.0],[15.588457268119894, 24.0],[20.784609690826528, 42.0],[20.784609690826528, 45.0],[18.186533479473212, 46.0],
#    [18.186533479473212, 49.0],[23.38268590217984, 46.0],[23.38268590217984, 49.0]]
#    K=V.coord[0:255]+V.coord[1536:1791] RECTANGULAR
    K=V.shell_cone
#    K=K.tolist()





for cell in Allcells.values(): 
     if cell.pos in K:
         cell.collagen_value=1


def find_key(dico,v):
    for k,val in dico.items():
        if v==val.cellid:
            return(k)
    return("la clef n'existe pas")
#collagen=[cell.pos for cell in Allcells.values() if cell.collagen_value==1]
#print(len(collagen),'len of collagen')
#for i in position_dico:
    #print(i)


#list_of_ids_collagen=[]
    

#for cell in Allcells.values():
    #if cell.pos in A:
        #list_of_ids_collagen.append(cell.cellid)


#for cell in Allcells.values():
 #   if cell.collagen_value==1:
  #      print('The cell has collagen')

#print('len of la liste',len(list_of_ids_collagen))



#for cell in Allcells.values():
    #if cell.cellid in list_of_ids_collagen:
        #cell.collagen_value=1



collagenpositions=[cell.pos for cell in Allcells.values() if cell.collagen_value==1]
for cell in Allcells.values():
    #if cell.collagen_value==1:
        collagenpositions.append(cell.pos)
        


liste_des_temps=[dt]
nombre_de_cellules=[len(V.coord)]
list_virus=[]
list_of_cancer_cell=[first_cancer_cell]
cancerous_cellbytime=[]
list_of_cancer_cell_infected=[]
Tmax=D.T_max
Already=False
cellule_deleted=[]
counter=0
while dt<Tmax:

    Positions_of_cells_newly_reproduced={}
    list_of_new_cancer_cells=[]
    w=[]
    z=[]
    for cell in Allcells.values():
        if cell.celltype==1 and cell.age(dt)>3:
            #print('La cellule est de type cancereuse et son age est >1, elle peut se reproduire')
            #print('La position de la cellule cancereuse est ',cell.pos)
            neighbors_list=cell.find_cellids_of_neighbors(V.k)
            #print('la liste des indices des voisins',neighbors_list)
            positions=[Allcells[i].pos for i in neighbors_list]
            distances=[cell.distance(Allcells[i]) for i in neighbors_list]
            #print('les differentes distances de la cellule par rapport aux voisins',distances)
            #print('la distance minimale autour de la cellule est ',min(distances))
            if len(positions)==0:
                #print('no neighbors available,it cannot reproduce anymore')
                continue
            if len(distances)!=0 and max(distances)>1.5 and min(distances)>0.5:
                #print('trop peu de place pour se reproduire')
                #print('la cellule peut se reproduire')
                #print('distances',distances)
                #print('la liste des indices des voisins',neighbors_list)
                #print('les differentes positions des voisins',positions)
                CA_collagen_values=[Allcells[i].collagen_value for i in neighbors_list]
                #print('CA_collagen_values',CA_collagen_values)
                f_CA_collagen_values=[g(k) for k in CA_collagen_values]
                #print('f_CA_collagen_values',f_CA_collagen_values)
                A=(10+sum(f_CA_collagen_values))/(1-np.exp(-rgrowth*dt))
                #A=10/1-np.exp(-rgrowth*dt)
                #print('A',A)
                CA_collagen_widths=[float(i/A) for i in f_CA_collagen_values]
                #print('CA collagen widths',CA_collagen_widths)
                CA_collagen_intervals=list_of_intervals(CA_collagen_widths)
                #print('CA_collagen_intervals',CA_collagen_intervals)
                p=random.random()
                if p<=CA_collagen_intervals[-1]:
                        #print("p est dans l'intervalle")
                        i=closest_element(CA_collagen_intervals,p)
                        #print('The index is',i)
                        neighbor_index=neighbors_list[i]
                        #print('It will follow the direction of this neighbor',neighbor_index)
                        neighbor=Allcells[neighbor_index]
                        #print('position du voisin suivi',neighbor.pos)
                        distance=cell.distance(Allcells[neighbor_index])
                        #print('la distance du voisin suivi est ' ,distance)
                        new_position=[cell.abscisse()+epsilon*cell.cos_theta(neighbor),cell.ordonnee()+epsilon*cell.sin_theta(neighbor),cell.profondeur()+epsilon*cell.Z_theta(neighbor)]
                        position_of_the_cell_who_reproduced=[cell.abscisse()-epsilon*cell.cos_theta(neighbor),cell.ordonnee()-epsilon*cell.sin_theta(neighbor),cell.profondeur()-epsilon*cell.Z_theta(neighbor)]
                        #print('Position de la cellule qui vient de naitre',new_position)
                        #print('nouvelle position de la cellule qui sest reproduite',position_of_the_cell_who_reproduced)
                        Positions_of_cells_newly_reproduced[cell.cellid-cell.shift]=position_of_the_cell_who_reproduced
                        new_cell=All_cells(getnext_cellid(),new_position,1,0,2.8,dt)
                        list_of_new_cancer_cells.append(new_cell)
                        cell.birthtime=dt
                        #print('la naissance de la cellule est ',cell.birthtime)
                        #print('lage de la cellule qui vient de se reproduire est',cell.age(dt))
                else:
                        #print("Il ne se passe rien, cette cellule ne se reproduit pas.")
                        continue
    bigkey=0
    for keys in Allcells.keys():
        if keys>bigkey:
            bigkey=keys
#    for kk in range(len(list_of_new_cancer_cells)):
#        list_of_new_cancer_cells[kk].shift=copy.deepcopy(Allcells[bigkey].shift)
    #ACTUALISATION DES SHIFTS
    for k in range (len(Allcells)):
        Allcells[k].shift=copy.deepcopy(Allcells[k].cellid-k)
    dico_nouvelles_positions={}

    for cell in Allcells.values():
        if cell.collagen_value==0:
            i=cell.cellid
            a,b,c=HC.nouvelle_position(i,V.k,1,Allcells)
            dico_nouvelles_positions[i-cell.shift]=[a,b,c]
        #print(i,E.coordonnees_bis[i],dico_nouvelles_positions[i])
        
    for i in dico_nouvelles_positions.keys():
        Allcells[i].pos=dico_nouvelles_positions[i]
        V.coord[i]=dico_nouvelles_positions[i]        
        
    for cell in list_of_new_cancer_cells:
        add_C_cell(cell)
        V.coord=np.append(V.coord,[list(cell.pos)],0)
      
    for id in Positions_of_cells_newly_reproduced.keys():
        #print('id de la cellule qui sest reproduite ',id )
        #print(Positions_of_cells_newly_reproduced[id])
        #print('E.coordonnees_bis[id]',E.coordonnees_bis[id])
        V.coord[id]=Positions_of_cells_newly_reproduced[id]
        #print('Changement effectue dans les coordonees des cellules qui viennent de se reproduire',E.coordonnees_bis[id])
        Allcells[id].pos=Positions_of_cells_newly_reproduced[id]
       # print('Changement effectue dans les coordonees des cellules qui viennent de se reproduire',Allcells[id].pos)

    collagenpositions=[cell.pos for cell in Allcells.values() if cell.collagen_value==1]
    #print(collagenpositions)
    list_of_cancer_cell=copy.deepcopy(list(set(list_of_new_cancer_cells+list_of_cancer_cell)))
    Xmin,Xmax=V.mini(V.coord)[0],V.maxi(V.coord)[0]
    Ymin,Ymax=V.mini(V.coord)[1],V.maxi(V.coord)[1]
    Zmin,Zmax=V.mini(V.coord)[2],V.maxi(V.coord)[2]
#    V.k=voroplusplus.compute_voronoi(V.coord,[[-40, 30.0],[-40, 30], [-15, 25]],1024.0,radii=[1.3, 1.4])
    V.k=copy.deepcopy(voroplusplus.compute_voronoi(V.coord,[[Xmin-1,Xmax+1],[Ymin-1, Ymax+1], [Zmin-1, Zmax+1]],1024.0,radii=[1.3, 1.4]))

    c1=copy.deepcopy(V.k)

 
    data.append(c1)
    pos_cancer_cell=copy.deepcopy([i.pos for i in list_of_cancer_cell]) 
    cancerous_cellbytime.append(pos_cancer_cell)
    resultat=[]
    if dt==1:
        density=D.diffusion(0,1)[0]
        for fa in range(D.Nx):
            for fa1 in range(D.Ny):
                for fa2 in range(D.Nz):
                    resultat.append(density[fa][fa1][fa2])
        grid=D.diffusion(0,1)[1]
    else:
        density=D.diffusion(dt-1,dt,density)[0]
        for fa in range(D.Nx):
            for fa1 in range(D.Ny):
                for fa2 in range(D.Nz):
                    resultat.append(density[fa][fa1][fa2])
    if (len(pos_cancer_cell)>0 and len(list_of_cancer_cell))>0:
        for i in range(len(pos_cancer_cell)):
            if resultat[closest_3Dposition(pos_cancer_cell[i],grid)]>0.005:
                
                
                
                for k in range (len(Allcells)):
                    Allcells[k].shift=copy.deepcopy(Allcells[k].cellid-k)
                for cancer_cell in list_of_cancer_cell:
                    for allcell in Allcells.values():
                        if allcell.cellid==cancer_cell.cellid:
                            cancer_cell.shift=copy.deepcopy(allcell.shift)
                for viruss in list_virus:
                    for cancer_cell in list_of_cancer_cell:
                        if viruss.cellid==cancer_cell.cellid:
                            viruss.shift=copy.deepcopy(cancer_cell.shift)
                            
                            
                            
                            
                for k in range (len(Allcells)):
                    Allcells[k].shift=copy.deepcopy(Allcells[k].cellid-k)
                nb_virus=resultat[closest_3Dposition(pos_cancer_cell[i],grid)]*c1[list_of_cancer_cell[i].cellid-list_of_cancer_cell[i].shift-1]['volume']
                if resultat[closest_3Dposition(pos_cancer_cell[i],grid)]>nb_virus/(D.hx*D.hy*D.hz):
                    resultat[closest_3Dposition(pos_cancer_cell[i],grid)]-=nb_virus/(D.hx*D.hy*D.hz)
                else:
                    resultat[closest_3Dposition(pos_cancer_cell[i],grid)]=0
                    il=0
                    for fa in range(D.Nx):
                        for fa1 in range(D.Ny):
                            for fa2 in range(D.Nz):
                                density[fa][fa1][fa2]=resultat[il]
                                il+=1
                Already=False
                for c in list_of_cancer_cell_infected:
                    if c.cellid==list_of_cancer_cell[i].cellid:
                        Already=True
 
                if Already==False:
                    list_virus.append(virus(list_of_cancer_cell[i].cellid,nb_virus,list_of_cancer_cell[i].shift))
                    list_virus[-1].birthtime=dt
                    list_of_cancer_cell_infected.append(list_of_cancer_cell[i])
                    print ("infection cancerous_cell",list_of_cancer_cell[i].cellid)
    

    iteration=0
        
    counter=0
    if len(list_virus)!=0:
        while iteration<len(list_virus):
            

#            for cancer_cell in list_of_cancer_cell:
#                for allcell in Allcells.values():
#                    if allcell.cellid==cancer_cell.cellid:
#                        cancer_cell.shift=copy.deepcopy(allcell.shift)
#            for viruss in list_virus:
#                for cancer_cell in list_of_cancer_cell:
#                    if viruss.cellid==cancer_cell.cellid:
#                        viruss.shift=copy.deepcopy(cancer_cell.shift)
            for k in range (len(Allcells)):
                Allcells[k].shift=copy.deepcopy(Allcells[k].cellid-k)
            idc=copy.deepcopy(list_virus[iteration].cellid)
            KEY=find_key(Allcells,idc)
            print(KEY)
            if list_virus[iteration].nb>10:
                if list_virus[iteration].birthtime<dt:
#                    if newindex==len(Allcells):
#                        cellule_deleted.append(Allcells[newindex])
#                        counterr+=1
#                    else:
#                        cellule_deleted.append(Allcells[newindex])
                    resultat[closest_3Dposition(Allcells[KEY].pos,grid)]=100
                    print("lysis cell",Allcells[KEY].cellid)
                    
                    
                    
                    
                    
                    del Allcells[KEY]
                    del list_virus[iteration]
                    iteration2=0
                    while iteration2 <len(list_of_cancer_cell_infected):
                        if list_of_cancer_cell_infected[iteration2].cellid==idc:
                            del list_of_cancer_cell_infected[iteration2]
                        iteration2+=1
                    iteration3=0
                    while  iteration3 <len(list_of_cancer_cell):
                            if list_of_cancer_cell[iteration3].cellid==idc:
                                del list_of_cancer_cell[iteration3]
                            iteration3+=1
                    print('the deleted key is:',KEY)
#                    
#                    pos_cancer_cell=[]
#                    for cell in list_of_cancer_cell:
#                        pos_cancer_cell.append(cell.pos)
#                    iteration4=0
#                    while  iteration4 <len(pos_cancer_cell):
#                        if pos_cancer_cell[iteration4]==celliddeleted:
#                            del pos_cancer_cell[iteration4]
#                        iteration4+=1
 
                    
                    
                    
                    for i in range(KEY,len(Allcells)):
                        Allcells[i]=Allcells.pop(i+1)
#                    iterationvirus=0
#                    while iterationvirus<len(list_virus):
#                        if list_virus[iterationvirus].cellid==cd:
#                            del list_virus[iterationvirus]
#                        iterationvirus+=1
#                    print("lysis cell",list_virus[iteration].cellid)
#                    
#                    
#                    
#                    
#                    iteration2=0





#                
                    il=0
                    for fa in range(D.Nx):
                        for fa1 in range(D.Ny):
                            for fa2 in range(D.Nz):
                                density[fa][fa1][fa2]=resultat[il]
                                il+=1

#                u=0
#                while u <len(list_of_cancer_cell):
#                    if list_of_cancer_cell[u].cellid==list_virus[iteration].cellid:
#                        del list_of_cancer_cell[u]
#                    u+=1
#                bu=0
                    new_coord=[]
                    for k in range(len(Allcells)):
                        new_coord.append(Allcells[k].pos)
                    V.coord=new_coord
#                for i in range(len(V.k)):
#                    for j in range(len(V.k[i]['faces'])):
#                        if V.k[i]['faces'][j]['adjacent_cell']>0:
#                            V.k[i]['faces'][j]['adjacent_cell']-=
                    V.k=voroplusplus.compute_voronoi(V.coord,[[Xmin-1,Xmax+1],[Ymin-1, Ymax+1], [Zmin-1, Zmax+1]],1024.0,radii=[1.3, 1.4])
                    test=V.k


            iteration+=1
    for viruss in list_virus : viruss.nb*=2
    if len(list_of_cancer_cell)==0:
        print("Cancer eliminated")
    dt+=1
    s+=1
    print('dt',dt-1)
    liste_des_temps.append(dt)
    nombre_de_cellules.append(len(V.coord))
#del Allcells[len(Allcells)-1]
data_bis=[collagenpositions,cancerous_cellbytime]
data=[data,data_bis]
fig,ax=plt.subplots()
ax.plot(liste_des_temps,nombre_de_cellules)
ax.set(xlabel='time (s)', ylabel='number of cells',title='Cell evolution with virus b=1000 gr=2 nb=2 mu=2')
ax.grid()

fig.savefig("cell(time)_with_virus_b=1000_gr=2__nb=2 mu=2.png")
plt.show()
print('liste des temps',liste_des_temps)
print('nombre de cellules',nombre_de_cellules)
#with open('data_cone_35.py','wb') as data1:
#    pickle.dump(data,data1)
#images=[f for f in os.listdir('.') if os.path.isfile(os.path.join('.'))]

