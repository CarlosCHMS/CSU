

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from su2MeshReader import reader
import sys
import numpy
import characteristics as ch
 
def csv2dict(fileName):
    data = numpy.genfromtxt(fileName, delimiter=',', names=True, dtype=None, encoding=None)
    
    return {name: data[name] for name in data.dtype.names}

class solution():

    def __init__(self, meshFile, solFile):
        
        self.mesh = reader(meshFile)
        
        self.x = self.mesh.x
        self.y = self.mesh.y
        self.elemToTri()
                
        self.pConnect()        
                
        ff = open(solFile)

        self.data = csv2dict(solFile)

        ff.close()
        
    def elemToTri(self):
    
        self.elem = []
        for e in self.mesh.elem:
            if(len(e) == 3):
                self.elem.append(e)
            elif(len(e) == 4):
                self.elem.append([e[0], e[1], e[2]])
                self.elem.append([e[2], e[3], e[0]])
            
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

    s = solution(path+"mesh.su2", path+"solution2.csv")

    triang = mtri.Triangulation(s.x, s.y, s.elem)

    plt.figure()
    plt.title("Static pressure")
    plt.tricontourf(triang, s.data['p'])
    #plt.triplot(triang, 'ko-') 
    plt.axis('equal') 
    plt.colorbar()  
    plt.show()

    plt.figure()
    plt.title("Mach")
    plt.tricontourf(triang, s.data['mach'])
    #plt.triplot(triang, 'ko-') 
    plt.axis('equal') 
    plt.colorbar()  
    plt.show()
        
    mar = s.mesh.markers[0]
    mar.getXY(s.mesh)
    
    inter = mtri.LinearTriInterpolator(triang, s.data['mach'])
    mach = inter(mar.x, mar.y)
    inter = mtri.LinearTriInterpolator(triang, s.data['p'])
    p = inter(mar.x, mar.y)
    inter = mtri.LinearTriInterpolator(triang, s.data['r'])
    r = inter(mar.x, mar.y)


    mar.x = numpy.array(mar.x)

    char1 = ch.problem(xchange=25.0, p1=1e4, T1=240, p4=1e5, T4=300)
    charVar = char1.calcVar(0.02, mar.x)

    plt.figure()
    plt.title("Mach")
    plt.plot(mar.x, mach, 'b')
    plt.plot(charVar['x'], charVar['u']/numpy.sqrt(1.4*charVar['p']/charVar['rho']),'--r')
    plt.show()
    
    plt.figure()
    plt.title("Static pressure")        
    plt.plot(mar.x, p, 'b')
    plt.plot(charVar['x'], charVar['p'],'--r')
    plt.show()
    
    plt.figure()
    plt.title("Density")        
    plt.plot(mar.x, r, 'b')
    plt.plot(charVar['x'], charVar['rho'],'--r')
    plt.show()
