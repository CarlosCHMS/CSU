

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from su2MeshReader import reader
import sys
import numpy
import readFluent as rf 
 
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
            
def qdin(p, T, m):

    gamma = 1.4
    Rgas = 287.0530
    
    r = p/(Rgas*T)
    c = numpy.sqrt(gamma*p/r)
    U = c*m
    
    return 0.5*r*U**2

    
            
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
    plt.tricontourf(triang, s.data['u'], levels(s.data['u'], 30))
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
    
    inter = mtri.LinearTriInterpolator(triang, s.data['u'])
    u = inter(mar.x, mar.y)
    inter = mtri.LinearTriInterpolator(triang, s.data['n'])
    n = inter(mar.x, mar.y)
    inter = mtri.LinearTriInterpolator(triang, s.data['T'])
    T = inter(mar.x, mar.y)
    inter = mtri.LinearTriInterpolator(triang, s.data['r'])
    r = inter(mar.x, mar.y)

    mar.y = numpy.array(mar.y)

    fm = rf.read(path+"resultsFluent/on_Ux")

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
    plt.plot(mar.y, T, '.b') 
    plt.show()

    fm = rf.read(path+"resultsFluent/on_shearStress")
    
    surf = csv2dict(path+'surfData.csv')

    plt.figure()
    plt.plot(surf['x'], -surf['Cfx'], 'b-') 
    plt.plot(fm.x, fm.y/qdin(1e5, 300, 0.1), 'r--')
    plt.legend(['Code', 'Fluent'])
    plt.xlabel("x [m/s]")
    plt.ylabel("Cfx [-]")
    plt.grid(True)
    plt.show()
    
    plt.figure()
    plt.plot(surf['x'], surf['yplus'], 'b-') 
    plt.xlabel("x [m/s]")
    plt.ylabel("yplus [-]")
    plt.grid(True)
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
