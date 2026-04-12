import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy
import sys
import os
import readFluent as rf 

        
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

    sys.path.append(os.path.abspath(path+'..'))
    from readers import CSUread

    s = CSUread(path)

    s.plotResiduals()
    
    s.plotConvergence('Cx_v')
    
    s.plotConvergence('Cy_p')

    s.plotSolution('p')
    
    s.plotSurfData('x', 'yplus')
    
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

    fm = rf.read(path+"resultsFluent/on_Ux")

    plt.figure()
    plt.plot(y, u, 'r--')    
    plt.plot(fm.x, fm.y, 'b')
    plt.legend(['Code', 'Fluent'])
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
    
    fm = rf.read(path+"resultsFluent/on_shearStress")
    
    plt.figure()
    plt.plot(s.surfData['x'], -s.surfData['Cfx'], 'b-') 
    plt.plot(fm.x, fm.y/qdin(1e5, 300, 0.1), 'r--')
    plt.legend(['Code', 'Fluent'])
    plt.xlabel("x [m/s]")
    plt.ylabel("Cfx [-]")
    plt.grid(True)
    plt.show()

    
