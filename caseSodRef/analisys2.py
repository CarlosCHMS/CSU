

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from su2MeshReader import reader
import sys
import numpy
import characteristics as ch
 
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
        
        for ii in range(0, len(self.elem)):
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

    s = solution(path+"mesh.su2", path+"solution1.csv")

    triang = mtri.Triangulation(s.x, s.y, s.elem)
        
    mar = s.mesh.markers[0]
    mar.getXY(s.mesh)
    
    inter = mtri.LinearTriInterpolator(triang, s.mach)
    mach1 = inter(mar.x, mar.y)
    inter = mtri.LinearTriInterpolator(triang, s.p)
    p1 = inter(mar.x, mar.y)
    inter = mtri.LinearTriInterpolator(triang, s.r)
    r1 = inter(mar.x, mar.y)

    s = solution(path+"mesh.su2", path+"solution2.csv")

    triang = mtri.Triangulation(s.x, s.y, s.elem)
        
    mar = s.mesh.markers[0]
    mar.getXY(s.mesh)
    
    inter = mtri.LinearTriInterpolator(triang, s.mach)
    mach2 = inter(mar.x, mar.y)
    inter = mtri.LinearTriInterpolator(triang, s.p)
    p2 = inter(mar.x, mar.y)
    inter = mtri.LinearTriInterpolator(triang, s.r)
    r2 = inter(mar.x, mar.y)

    mar.x = numpy.array(mar.x)

    char1 = ch.problem(xchange=25.0, p1=1e4, T1=240, p4=1e5, T4=300)
    charVar = char1.calcVar(0.02, mar.x)

    plt.figure()
    plt.plot(mar.x, mach1, 'b')
    plt.plot(mar.x, mach2, 'g')
    plt.plot(charVar['x'], charVar['u']/numpy.sqrt(1.4*charVar['p']/charVar['rho']),'--r')
    plt.legend(['1 order', '2 order', 'charac.'])
    plt.grid(True)
    plt.xlabel("x [m]")
    plt.ylabel("mach [-]")
    plt.savefig('sod_u.png', dpi=300)
    
    plt.figure()
    plt.title("Static pressure")        
    plt.plot(mar.x, p1, 'b')
    plt.plot(mar.x, p2, 'g')
    plt.plot(charVar['x'], charVar['p'],'--r')
    plt.legend(['1 order', '2 order', 'charac.'])
    plt.grid(True)
    plt.savefig('sod_pressure.png', dpi=300)
    
    plt.figure()      
    plt.plot(mar.x, r1, 'b')
    plt.plot(mar.x, r2, 'g')
    plt.plot(charVar['x'], charVar['rho'],'--r')
    plt.legend(['1 order', '2 order', 'charac.'])
    plt.grid(True)
    plt.xlabel("x [m]")
    plt.ylabel(r"$\rho$ [kg/m³]")    
    plt.savefig('sod_rho.png', dpi=300)

