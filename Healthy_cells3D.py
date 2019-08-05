import numpy as np 
import matplotlib.pyplot as plt 
import voronoi3D as V
import scipy.spatial as sp 
import os 
import random

H_C={}

class Healthy_cell():
    
    def __init__(self,cellid,pos,celltype=0,collagen_value=0,spring=2.8):
        self.cellid=cellid
        self.pos=pos
        self.celltype=celltype
        self.collagen_value=collagen_value
        self.spring=spring

    def distance(self,other):
        return np.sqrt((self.pos[0]-other.pos[0])**2+(self.pos[1]-other.pos[1])**2+(self.pos[2]-other.pos[2])**2)
    
    def distance_x(self,other):
        return(self.pos[0]-other.pos[0])
    
    def distance_y(self,other):
        return(self.pos[1]-other.pos[1])
        
    def distance_z(self,other):
        return(self.pos[2]-other.pos[2])

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
    
    def find_cellids_of_neighbors(self,k):
        #returns a list of all the cell_ids of the neighbors.
        neighbors=[]
        point_index=self.cellid
        cell=k[point_index]
        for faces in cell["faces"]:
            neighbors.append(faces['adjacent_cell'])
        return list(set(neighbors))




for i in range(len(V.coord)):
    H_C[i]=Healthy_cell(i,(V.coord[i][0],V.coord[i][1],V.coord[i][2]))



#for i in E.coordonnees_bis:
    #print(i)

#for i in range(len(H_C)):
    #cell=H_C[i]
    #print(cell.cellid,cell.pos,cell.find_cellids_of_neighbors(E.tri))

def nouvelle_position(i,k,dt,dico):
    eta=1
    for k1 in range (len(dico)):
        if dico[k1].cellid==i:
            our_cell=dico[k1]
    neighbor_list_filtree=[]
    Q=0
    neighbor_list=our_cell.find_cellids_of_neighbors(k)
    xm,ym,zm=0,0,0
    distances=[our_cell.distance(dico[j]) for j in neighbor_list]
    for i in range(len(neighbor_list)):
        if distances[i]<3:
            neighbor_list_filtree.append(neighbor_list[i])
    #neighbor_list_filtree=[j for j in neighbor_list if distances[j]<3]
    if our_cell.celltype==0:
        mu=0.1
        x=0
        y=0
        z=0
        for j in neighbor_list_filtree:
            neighbor_cell=dico[j]
            distance=our_cell.distance(neighbor_cell)
            coeff_cancer_healthy=1

            if neighbor_cell.celltype==1:
                 coeff_cancer_healthy=5
            if our_cell.abscisse()-neighbor_cell.abscisse()<0:
                x+=-mu*(neighbor_cell.spring-distance)*coeff_cancer_healthy
            else:
                x+=mu*(neighbor_cell.spring-distance)*coeff_cancer_healthy
            if our_cell.ordonnee()-neighbor_cell.ordonnee()>0:
                y+=mu*(neighbor_cell.spring-distance)*coeff_cancer_healthy
            else:
                y+=-mu*(neighbor_cell.spring-distance)*coeff_cancer_healthy
            if our_cell.profondeur()-neighbor_cell.profondeur()>0:
                z+=mu*(neighbor_cell.spring-distance)*coeff_cancer_healthy
            else:
                z+=-mu*(neighbor_cell.spring-distance)*coeff_cancer_healthy

            
        F_x=x
        F_y=x
        F_z=z
        new_x=our_cell.abscisse()+dt*F_x/eta

        new_y=our_cell.ordonnee()+dt*F_y/eta

        new_z=our_cell.profondeur()+dt*F_z/eta
         
    if our_cell.celltype==1:
        mu=0.1
        x,y,z=0,0,0
        new_x,new_y,new_z=our_cell.abscisse(),our_cell.ordonnee(),our_cell.profondeur()
        if len(neighbor_list_filtree)!=0:
            for j in neighbor_list_filtree:
                coeff_cancer_healthy=1
                coeff_cancer_collagen=1
                neighbor_cell=dico[j]
                distance=our_cell.distance(neighbor_cell)
                if neighbor_cell.celltype==0:
                        coeff_cancer_healthy=5
                if neighbor_cell.collagen_value==1:
                        coeff_cancer_collagen=2
                #print(distance)
                #print(neighbor_cell.spring))
                if our_cell.abscisse()-neighbor_cell.abscisse()<0:
                    x+=-mu*(neighbor_cell.spring-distance)*coeff_cancer_healthy*coeff_cancer_collagen
                else:
                    x+=mu*(neighbor_cell.spring-distance)*coeff_cancer_healthy*coeff_cancer_collagen
                if our_cell.ordonnee()-neighbor_cell.ordonnee()>0:
                    y+=mu*(neighbor_cell.spring-distance)*coeff_cancer_healthy*coeff_cancer_collagen
                else:
                    y+=-mu*(neighbor_cell.spring-distance)*coeff_cancer_healthy*coeff_cancer_collagen
                if our_cell.profondeur()-neighbor_cell.profondeur()>0:
                    z+=mu*(neighbor_cell.spring-distance)*coeff_cancer_healthy*coeff_cancer_collagen
                else:
                    z+=-mu*(neighbor_cell.spring-distance)*coeff_cancer_healthy*coeff_cancer_collagen
                ###### METASTASIS ATTEMPS
#                xm,ym,zm=x,y,z
            
#            if randomy>0.995 :
#                x,y=0,0
#                F_x=x+2*xm*random.random()
#                #print('F_x:',F_x)
#                F_y=y+2*ym*random.random()
#                F_z=z+2*zm*random.random()
#                our_cell.celltype=2
#                Q=dt
#            else:
#                F_x=x
#                F_y=y
#                F_z=z
            #########

            F_x=x
            F_y=y
            F_z=z
            new_x=our_cell.abscisse()+dt*F_x/eta
            new_y=our_cell.ordonnee()+dt*F_y/eta
            new_z=our_cell.profondeur()+dt*F_z/eta
            
    #####METASTASIS_ATTEMP
#    if our_cell.celltype==2:
#        new_x=our_cell.abscisse()+dt*1*random.random()
#        new_y=our_cell.ordonnee()+dt*1*random.random()
#        new_z=our_cell.profondeur()+dt*1*random.random()
#        if dt>Q+10:
#            our_cell.celltype=1
    ########
    return new_x,new_y,new_z


