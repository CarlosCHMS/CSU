

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from su2MeshReader import reader
import sys
import numpy
import readFluent as rf 
 
class solution():

    def __init__(self, meshFile, solFile):

        self.gamma = 1.4
        self.Rgas = 287.5

        self.mesh = reader(meshFile)
        
        self.x = self.mesh.x
        self.y = self.mesh.y
        self.elemToTri()
                
        self.pConnect()        
                
        ff = open(solFile)

        self.r = []
        self.ru = []
        self.rv = []
        self.rE = []
        self.rn = []

        first = True
        for row in ff:
            if first:
                first = False
            else:
                aux = row.split(',')
                self.r.append(float(aux[0]))
                self.ru.append(float(aux[1]))
                self.rv.append(float(aux[2]))
                self.rE.append(float(aux[3]))
                if len(aux) == 6:
                    self.flagN = True
                    self.rn.append(float(aux[4]))

        ff.close()

        self._toArray()
        self.calcPMT()
    
    def _toArray(self):

        self.r = numpy.array(self.r)
        self.ru = numpy.array(self.ru)
        self.rv = numpy.array(self.rv)
        self.rE = numpy.array(self.rE)

        return None
        
    def elemToTri(self):
    
        self.elem = []
        for e in self.mesh.elem:
            if(len(e) == 3):
                self.elem.append(e)
            elif(len(e) == 4):
                self.elem.append([e[0], e[1], e[2]])
                self.elem.append([e[2], e[3], e[0]])
            
        return None                
        
    def calcPMT(self):
       
        Np = len(self.mesh.p)
       
        self.p = numpy.zeros(Np)
        self.mach = numpy.zeros(Np)
        self.entro = numpy.zeros(Np)               
        self.H = numpy.zeros(Np)                
        self.u = numpy.zeros(Np)
        self.T = numpy.zeros(Np)
        self.n = numpy.zeros(Np)
       
        for ii in range(0, Np):
       
            if(self.con[ii] > 0):
       
                u = self.ru[ii]/self.r[ii]
                v = self.rv[ii]/self.r[ii]
                E = self.rE[ii]/self.r[ii]
                self.n[ii] = self.rn[ii]/self.r[ii]                

                self.u[ii] = u
                
                RT = (E - (u**2 + v**2)/2)*(self.gamma - 1)
                
                self.p[ii] = RT*self.r[ii]
                
                c = numpy.sqrt(self.gamma*RT)
                V = numpy.sqrt(u**2 + v**2)
                
                self.mach[ii] = V/c

                self.entro[ii] = self.p[ii]/(self.r[ii]**self.gamma)

                self.H[ii] = E + self.p[ii]/self.r[ii]
                
                self.T[ii] = RT/self.Rgas
        
        return None        

    def pConnect(self):
    
        self.con = numpy.zeros(len(self.mesh.p))
        
        for ii in range(0, len(self.elem)):
            self.con[self.elem[ii][0]] += 1
            self.con[self.elem[ii][1]] += 1
            self.con[self.elem[ii][2]] += 1
            
        return None

def levels(v, n):    

    max1 = v[0]
    min1 = v[0]
    for ii in range(0, v.shape[0]):
        max1 = max(v[ii], max1)
        min1 = min(v[ii], min1)
                            
    d = (max1-min1)/(n-1)
    levels = []
    for ii in range(0, n):
        levels.append(min1 + d*ii)
    
    return levels                
    
    
class BL():

    def __init__(self, p, T, m, L):
    
        # Blasius solution
    
        gamma = 1.4
        Rgas = 287.5
        
        r = p/(Rgas*T)
        self.c = numpy.sqrt(gamma*p/r)
        self.U = self.c*m
        
        mi = 1.45e-6*T*numpy.sqrt(T)/(T + 110.0)
        
        Re = r*self.U*L/mi
        
        self.h = L/numpy.sqrt(Re)

        self.tab = [[0, 0],
                    [0.5, 0.16503],
                    [1, 0.32819], 
                    [1.5, 0.48471],
                    [2, 0.62755 ],
                    [2.5, 0.74927],
                    [3, 0.84452 ],
                    [3.5,  0.91205],
                    [4, 0.95499 ],
                    [4.5, 0.97929 ],
                    [4.91, 0.98991 ],
                    [4.92, 0.99009 ],
                    [5, 0.99147 ],
                    [6, 0.99898 ],
                    [7, 0.99993 ],
                    [8, 1]]

        for t in self.tab:
            t = numpy.array(t)
            
        self.tab = numpy.array(self.tab)

        self.y = self.h*self.tab[:, 0]        
        self.u = self.U*self.tab[:, 1]                      
    
class convergence():

    def __init__(self, convFile):
    
        ff = open(convFile, 'r')
        ii = 0
        self.varList = []
        for row in ff:
            aux = row.split(',')
            if ii == 0:
                for jj in range(len(aux)-1):
                    self.varList.append([])

            for jj in range(len(aux)-1):
                self.varList[jj].append(float(aux[jj]))
                        
            ii += 1
    
        for jj in range(len(aux)-1):
            self.varList[jj] = numpy.array(self.varList[jj])  
            
            
if __name__=="__main__":

    if len(sys.argv) < 1:
        print("Insert the input file")
        exit(0)
    
    path = sys.argv[1]



    s = solution(path+"mesh.su2", path+"solution.csv")

    triang = mtri.Triangulation(s.x, s.y, s.elem)

    plt.figure()
    plt.title("Static pressure")
    plt.tricontourf(triang, s.p)
    #plt.triplot(triang, 'ko-') 
    plt.axis('equal') 
    plt.colorbar()  
    plt.show()

    plt.figure()
    #plt.title("Mach")
    plt.tricontourf(triang, s.u, levels=levels(s.u, 30))
    #plt.triplot(triang, 'ko-') 
    #plt.axis('equal') 
    cbar = plt.colorbar()
    cbar.set_label('u [m/s]')
    plt.xlabel("x [m]")  
    plt.ylabel("y [m]")    
    #plt.savefig("u.png", dpi=300)
    plt.show()
        
    mar = s.mesh.markers[1]
    mar.getXY(s.mesh)
    
    inter = mtri.LinearTriInterpolator(triang, s.u)
    u = inter(mar.x, mar.y)
    inter = mtri.LinearTriInterpolator(triang, s.n)
    n = inter(mar.x, mar.y)
    inter = mtri.LinearTriInterpolator(triang, s.T)
    T = inter(mar.x, mar.y)
    inter = mtri.LinearTriInterpolator(triang, s.r)
    r = inter(mar.x, mar.y)

    bl = BL(1e5, 300, 0.1, 0.5)

    mar.y = numpy.array(mar.y)

    u[0] = 0;


    fm = rf.read(path+"BL_turb/on_Ux")

    plt.figure()
    plt.plot(mar.y, u, 'b')
#    plt.plot(bl.y, bl.u, 'r--')
    plt.plot(fm.x, fm.y, 'g.')   
    plt.legend(['Code', 'Fluent'])
    plt.grid(True)    
    plt.xlabel("y [m]")  
    plt.ylabel("u [m/s]")        
#    plt.savefig(path+"uProfile.png", dpi=300)
    plt.show()    


    plt.figure()
    plt.title("T")
    plt.plot(mar.y, T, 'b') 
    plt.show()

    fm = rf.read(path+"BL_turb/on_ni")

    plt.figure()
    plt.plot(mar.y, n, 'b')
#    plt.plot(bl.y, bl.u, 'r--')
    plt.plot(fm.x, fm.y, 'g.')   
    plt.legend(['Code', 'Fluent'])
    plt.grid(True)    
    plt.xlabel("y [m]")  
    plt.ylabel("ni")        
    plt.show()

    
    conv = convergence(path+"convergence.csv")
    
    plt.figure()
    plt.semilogy(conv.varList[1]/conv.varList[1][0])
    plt.semilogy(conv.varList[2]/conv.varList[2][0])    
    plt.semilogy(conv.varList[3]/conv.varList[3][0])    
    plt.semilogy(conv.varList[4]/conv.varList[4][0])       
    plt.semilogy(conv.varList[5]/conv.varList[5][0])           
    plt.grid(True)
    plt.xlabel("iterations [100]")  
    plt.ylabel("Residuos")            
    plt.show()
    
    plt.figure()
    plt.plot(conv.varList[7])
    plt.grid(True)
    plt.xlabel("iterations [100]")  
    plt.ylabel("Cx_v [-]")            
    plt.show()
    
    plt.figure()
    plt.plot(conv.varList[8])
    plt.grid(True)
    plt.xlabel("iterations [100]")  
    plt.ylabel("Cy_p [-]")            
    plt.show()    
