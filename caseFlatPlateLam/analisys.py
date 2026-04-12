

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
    
    
class BL():

    def __init__(self, p, T, m, L):
    
        # Blasius solution
    
        gamma = 1.4
        Rgas = 287.0530
        
        r = p/(Rgas*T)
        self.c = numpy.sqrt(gamma*p/r)
        self.U = self.c*m
        
        self.qdin = 0.5*r*self.U**2
        
        mi = 1.45e-6*T*numpy.sqrt(T)/(T + 110.0)

        self.aux = 0.332*numpy.sqrt(r*mi*self.U**3)/self.qdin
        
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
        
        return None

    def calcCfx(self, x):

        return self.aux/numpy.sqrt(x)
        
    
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
    #plt.title("Mach")
    plt.tricontourf(triang, s.data['mach'])
    #plt.triplot(triang, 'ko-') 
    #plt.axis('equal') 
    plt.colorbar()  
    plt.show()

    plt.figure()
    #plt.title("Mach")
    plt.tricontourf(triang, s.data['u'])
    #plt.triplot(triang, 'ko-') 
    #plt.axis('equal') 
    cbar = plt.colorbar()
    cbar.set_label('u [m/s]')
    plt.xlabel("x [m]")  
    plt.ylabel("y [m]")
    plt.show()

    mar = s.mesh.markers[1]
    mar.getXY(s.mesh)
    
    inter = mtri.LinearTriInterpolator(triang, s.data['u'])
    u = inter(mar.x, mar.y)
    inter = mtri.LinearTriInterpolator(triang, s.data['T'])
    T = inter(mar.x, mar.y)
    inter = mtri.LinearTriInterpolator(triang, s.data['r'])
    r = inter(mar.x, mar.y)

    y = numpy.array(mar.y)

    bl = BL(1e5, 300, 0.1, 0.5)

    plt.figure()
    plt.plot(y, u, 'r--')    
    plt.plot(bl.y, bl.u, 'b')
    plt.legend(['Code', 'Blasius'])
    plt.grid(True)    
    plt.xlabel("y [m]")  
    plt.ylabel("u [m/s]")
    plt.show()
    
    plt.figure()
    plt.title("T")
    plt.plot(y, T, 'b') 
    plt.show()

    surf = csv2dict(path+'surfData.csv')

    plt.figure()
    plt.plot(surf['x'], -surf['Cfx'], 'b-') 
    plt.plot(surf['x'], bl.calcCfx(surf['x']), 'r--')
    plt.legend(['Code', 'Blasius'])    
    #plt.ylim([0, 7])    
    plt.xlabel("x [m/s]")
    plt.ylabel("Cfx [-]")
    plt.grid(True)
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
    plt.plot(conv['Cx_v'])
    plt.grid(True)
    plt.xlabel("iterations")  
    plt.ylabel("Cx_v")            
    plt.show()
