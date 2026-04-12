import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import sys
import os
import numpy
import characteristics as ch
        

if __name__=="__main__":

    if len(sys.argv) < 1:
        print("Insert the input file")
        exit(0)
    
    path = sys.argv[1]

    sys.path.append(os.path.abspath(path+'..'))
    from readers import CSUread

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
