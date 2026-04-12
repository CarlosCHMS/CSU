import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import sys
import os
import numpy
import characteristics as ch

class CSUread():

    def __init__(self, path):
        
        self.mesh = reader(path+'mesh.su2')

        self.solution = self.csv2dict(path+'solution2.csv')
        
        self.surfData = self.csv2dict(path+'surfData.csv')
        
        self.convergence = self.csv2dict(path+'convergence.csv')
        
        self.mesh.removeExtraPoints()

        if self.mesh.extraPoints:
            for k in self.solution.keys():
                self.solution[k] = self.solution[k][self.mesh.newNodes]        
        
        self.elemToTri()
        
        self.triang = mtri.Triangulation(self.mesh.x, self.mesh.y, self.elem)
        
    def elemToTri(self):
    
        self.elem = []
        for e in self.mesh.elem:
            if(len(e) == 3):
                self.elem.append(e)
            elif(len(e) == 4):
                self.elem.append([e[0], e[1], e[2]])
                self.elem.append([e[2], e[3], e[0]])
            
        return None
        
    def csv2dict(self, fileName):
        
        data = numpy.genfromtxt(fileName, delimiter=',', names=True, dtype=None, encoding=None)
        
        return {name: data[name] for name in data.dtype.names}                
        
    def plotSolution(self, field):
    
        plt.figure()
        plt.title(field)
        plt.tricontourf(self.triang, self.solution[field], levels=30)
        #plt.triplot(triang, 'ko-') 
        plt.axis('equal') 
        cbar = plt.colorbar()  
        cbar.set_label(field)
        plt.xlabel('x')
        plt.ylabel('y')        
        plt.show()

        return None
        
    def plotResiduals(self):
    
        plt.figure()
        leg = []        
        for k in s.convergence.keys():
            if 'res_' in k:
                plt.semilogy(s.convergence[k]/s.convergence[k][0])
                leg.append(k)

        plt.grid(True)
        plt.xlabel("iterations")  
        plt.ylabel("residuals")            
        plt.legend(leg)
        plt.show()
        
        return None
        
    def plotConvergence(self, field):
    
        plt.figure()
        plt.plot(s.convergence[field])
        plt.grid(True)
        plt.xlabel("iterations")  
        plt.ylabel(field)
        plt.show()
        
        return None        

    def plotSurfData(self, field1, field2):
    
        plt.figure()
        plt.plot(s.surfData[field1], s.surfData[field2])
        plt.grid(True)
        plt.xlabel(field1)  
        plt.ylabel(field2)
        plt.show()
        
        return None
        

if __name__=="__main__":

    if len(sys.argv) < 1:
        print("Insert the input file")
        exit(0)
    
    path = sys.argv[1]

    sys.path.append(os.path.abspath(path+'..'))
    from su2MeshReader import reader

    s = CSUread(path)

    s.plotSolution('p')
    
    s.plotSolution('mach')
        
    mar = s.mesh.markers[0]
    mar.getXY(s.mesh)
    
    mar.x = numpy.array(mar.x)

    char1 = ch.problem(xchange=25.0, p1=1e4, T1=240, p4=1e5, T4=300)
    charVar = char1.calcVar(0.02, mar.x)

    plt.figure()
    plt.title("Mach")
    plt.plot(s.surfData['x'], s.surfData['mach'], 'b')
    plt.plot(charVar['x'], charVar['u']/numpy.sqrt(1.4*charVar['p']/charVar['rho']),'--r')
    plt.show()
    
    plt.figure()
    plt.title("Static pressure")        
    plt.plot(s.surfData['x'], s.surfData['p'], 'b')
    plt.plot(charVar['x'], charVar['p'],'--r')
    plt.show()
    
    plt.figure()
    plt.title("Density")        
    plt.plot(s.surfData['x'], s.surfData['rho'], 'b')
    plt.plot(charVar['x'], charVar['rho'],'--r')
    plt.show()
