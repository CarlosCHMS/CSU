

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from su2MeshReader import reader
import sys
import numpy
 
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

    
if __name__=="__main__":

    if len(sys.argv) < 1:
        print("Insert the input file")
        exit(0)
    
    path = sys.argv[1]

    s = solution(path+"mesh.su2", path+"solution2.csv")

    triang = mtri.Triangulation(s.x, s.y, s.elem)

    plt.figure()
    plt.title("Static pressure")
    plt.tricontourf(triang, s.data['p'], levels=levels(s.data['p'], 20))
    #plt.triplot(triang, 'ko-') 
    plt.axis('equal') 
    plt.colorbar()  
    plt.show()

    plt.figure()
    plt.title("Mach")
    plt.tricontourf(triang, s.data['mach'], levels=levels(s.data['mach'], 20))
    #plt.triplot(triang, 'ko-') 
    plt.axis('equal') 
    plt.colorbar()  
    plt.show()
        
    conv = csv2dict((path+"convergence.csv"))
    
    plt.figure()
    plt.semilogy(conv['res_r']/conv['res_r'][0])
    plt.semilogy(conv['res_u']/conv['res_u'][0])
    plt.semilogy(conv['res_v']/conv['res_v'][0])
    plt.semilogy(conv['res_E']/conv['res_E'][0])
    plt.grid(True)
    plt.xlabel("iterations")  
    plt.ylabel("residuals")            
    plt.show()
    
    plt.figure()
    plt.plot(conv['Cx_p'])
    plt.plot(conv['Cy_p'])    
    plt.grid(True)
    plt.xlabel("iterations")
    plt.ylabel("coeff")
    plt.legend(["Cx_p", "Cy_p"])
    plt.show()
        
