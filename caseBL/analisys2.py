

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from su2MeshReader import reader
import sys
import numpy
import readFluent as rf 
 
class friction():

    def __init__(self):
    
        ff = open("./friction.csv", 'r')
        self.x = []
        self.txy = []
        self.yp = []
        for row in ff:
            aux = row.split(',')
            self.x.append(float(aux[0]))
            self.txy.append(float(aux[2]))
            self.yp.append(float(aux[4]))
            
        self.x = numpy.array(self.x)
        self.yp = numpy.array(self.yp)        
        
    
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

        self.aux = 0.332*numpy.sqrt(r*mi*self.U**3)
        
        self.qdin = 0.5*r*self.U**2

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
    
    def calcTxy(self, x):
    
        return self.aux/numpy.sqrt(x)
        
    def trapz(self, x, y):
    
        ans = 0.0
        for ii in range(0, len(x)-1):
            ans += 0.5*(x[ii+1] - x[ii])*(y[ii+1] + y[ii])
            
        return ans
            
    
if __name__=="__main__":

    f = friction()
    
    bl = BL(1e5, 300, 0.1, 0.5)
    
    print(bl.trapz(f.x, f.txy)/bl.qdin, bl.trapz(f.x, bl.calcTxy(f.x))/bl.qdin)    
    
    fm = rf.read("./frictionBLfluentRef")
    
    plt.figure()
    plt.plot(f.x, f.txy, 'b') 
    plt.plot(f.x, bl.calcTxy(f.x), 'r--')
    plt.legend(['Code', 'Blasius'])    
    plt.ylim([0, 7])    
    plt.xlabel("x [m/s]")
    plt.ylabel(r"$\tau_{xy}$ [Pa]")
    plt.grid(True)
    plt.savefig("./shearStress.png", dpi=300)
    
    plt.figure()
    plt.plot(f.x, f.yp,'--')
    plt.plot(f.x, 0.1+0*f.yp,'--')    
    plt.show()    
