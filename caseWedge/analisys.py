

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from su2MeshReader import reader
import sys
import numpy
 
class solution():

    def __init__(self, meshFile, solFile):

        self.gamma = 1.4
        self.Rgas = 287.5

        r = reader(meshFile)
        
        self.x = r.x
        self.y = r.y
        self.elem = r.elem
                
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
       
        u = self.ru/self.r
        v = self.rv/self.r
        E = self.rE/self.r
        
        RT = (E - (u**2 + v**2)/2)*(self.gamma - 1)
        self.T = RT/self.Rgas
        
        self.p = RT*self.r
        
        c = numpy.sqrt(self.gamma*RT)
        V = numpy.sqrt(u**2 + v**2)
        
        self.mach = V/c
        self.u = u        
        self.v = v

        self.entro = self.p/(self.r**self.gamma)

        self.H = E + self.p/self.r
        
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
