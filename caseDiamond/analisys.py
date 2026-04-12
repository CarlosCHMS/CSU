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
    
    s.plotSolution('s')
    
    s.plotSolution('H')

    s.plotSurfData('x', 'Cp')

