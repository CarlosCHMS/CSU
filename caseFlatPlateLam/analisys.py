import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy
import sys
import os
        
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

    sys.path.append(os.path.abspath(path+'..'))
    from readers import CSUread

    s = CSUread(path)

    s.plotResiduals()
    
    s.plotConvergence('Cx_v')
    
    s.plotConvergence('Cy_p')

    s.plotSolution('p')
    
    plt.figure()
    plt.title('mach')
    plt.tricontourf(s.triang, s.solution['mach'], levels=30)
    #plt.triplot(triang, 'ko-') 
    cbar = plt.colorbar()  
    cbar.set_label('mach')
    plt.xlabel('x')
    plt.ylabel('y')        
    plt.show()
    
    mar = s.mesh.markers[1]
    mar.getXY(s.mesh)
    
    inter = mtri.LinearTriInterpolator(s.triang, s.solution['u'])
    u = inter(mar.x, mar.y)
    inter = mtri.LinearTriInterpolator(s.triang, s.solution['T'])
    T = inter(mar.x, mar.y)

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
    plt.xlabel("y [m]")  
    plt.ylabel("T [K]")
    plt.show()    
    
    plt.figure()
    plt.plot(s.surfData['x'], -s.surfData['Cfx'], 'b-') 
    plt.plot(s.surfData['x'], bl.calcCfx(s.surfData['x']), 'r--')
    plt.legend(['Code', 'Blasius'])
    plt.xlabel("x [m/s]")
    plt.ylabel("Cfx [-]")
    plt.grid(True)
    plt.show()


