

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
        self.elem = self.mesh.elem
                
        self.pConnect()        
                
        ff = open(solFile)

        self.r = []
        self.ru = []
        self.rv = []
        self.rE = []

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

        ff.close()

        self._toArray()
        self.calcPMT()
    
    def _toArray(self):

        self.r = numpy.array(self.r)
        self.ru = numpy.array(self.ru)
        self.rv = numpy.array(self.rv)
        self.rE = numpy.array(self.rE)

        return None
        
    def calcPMT(self):
       
        Np = len(self.mesh.p)
       
        self.p = numpy.zeros(Np)
        self.mach = numpy.zeros(Np)
        self.entro = numpy.zeros(Np)               
        self.H = numpy.zeros(Np)                
       
        for ii in range(0, Np):
       
            if(self.con[ii] > 0):
       
                u = self.ru[ii]/self.r[ii]
                v = self.rv[ii]/self.r[ii]
                E = self.rE[ii]/self.r[ii]
                
                RT = (E - (u**2 + v**2)/2)*(self.gamma - 1)
                
                self.p[ii] = RT*self.r[ii]
                
                c = numpy.sqrt(self.gamma*RT)
                V = numpy.sqrt(u**2 + v**2)
                
                self.mach[ii] = V/c

                self.entro[ii] = self.p[ii]/(self.r[ii]**self.gamma)

                self.H[ii] = E + self.p[ii]/self.r[ii]
        
        return None        

    def pConnect(self):
    
        self.con = numpy.zeros(len(self.mesh.p))
        
        for ii in range(0, len(self.mesh.elem)):
            self.con[self.elem[ii][0]] += 1
            self.con[self.elem[ii][1]] += 1
            self.con[self.elem[ii][2]] += 1
            
        return None

def levels(v, n):    

    max1 = v[0][0]
    min1 = v[0][0]
    for ii in range(0, v.shape[0]):
        for jj in range(0, v.shape[1]):
            max1 = max(v[ii][jj], max1)
            min1 = min(v[ii][jj], min1)
                            
    d = (max1-min1)/(n-1)
    levels = []
    for ii in range(0, n):
        levels.append(min1 + d*ii)
    
    return levels                
    
if __name__=="__main__":

    if len(sys.argv) < 1:
        print("Insert the input file")
        exit(0)
    
    path = sys.argv[1]

    s = solution(path+"mesh.su2", path+"solution.csv")

    triang = mtri.Triangulation(s.x, s.y, s.elem)

    plt.figure()
    plt.tricontourf(triang, s.r)
    #plt.triplot(triang, 'ko-') 
    plt.axis('equal') 
    plt.colorbar()  
    plt.show()

    plt.figure()
    plt.tricontourf(triang, s.mach)
    #plt.triplot(triang, 'ko-') 
    plt.axis('equal') 
    plt.colorbar()  
    plt.show()
        
    mar = s.mesh.markers[0]
    mar.getXY(s.mesh)
    
    inter = mtri.LinearTriInterpolator(triang, s.mach)
    mach = inter(mar.x, mar.y)
    inter = mtri.LinearTriInterpolator(triang, s.p)
    p = inter(mar.x, mar.y)

    fm = rf.read(path+"machxX")

    mar.x = numpy.array(mar.x)

    plt.figure()
    plt.plot(mar.x, mach)
    plt.plot(mar.x, mar.x*0 + 2.0)    
    plt.plot(mar.x, mar.x*0 + 1.21021838)    
    plt.plot(fm.x, fm.y)
    plt.show()
    
    fp = rf.read(path+"pressurexX")    
    
    plt.figure()
    plt.plot(mar.x, p)
    plt.plot(mar.x, mar.x*0 + 1e5)    
    plt.plot(mar.x, mar.x*0 + 2.84286270*1e5)    
    plt.plot(fp.x, fp.y)    
    plt.show()
    
