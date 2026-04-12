import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import sys
import os
import numpy
        

if __name__=="__main__":

    if len(sys.argv) < 1:
        print("Insert the input file")
        exit(0)
    
    path = sys.argv[1]

    sys.path.append(os.path.abspath(path+'..'))
    from readers import CSUread

    s = CSUread(path)

    s.plotResiduals()
    
    s.plotConvergence('Cx_p')
    
    s.plotConvergence('Cy_p')

    s.plotSolution('p')
    
    s.plotSolution('mach')

    s.plotSurfData('x', 'Cp')
    
    x = s.surfData['x']
    
    plt.figure()
    plt.title("Mach")
    plt.plot(x, s.surfData['mach'])
    plt.plot(x, x*0 + 2.0)    
    plt.plot(x, x*0 + 1.21021838)    
    plt.show()
    
    plt.figure()
    plt.title("Static pressure")        
    plt.plot(x, s.surfData['p'])
    plt.plot(x, x*0 + 1e5)    
    plt.plot(x, x*0 + 2.84286270*1e5)    
    plt.show()    
    
    
