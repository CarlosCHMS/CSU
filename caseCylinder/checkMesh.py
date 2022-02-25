

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from su2MeshReader import reader
import sys
import numpy
               
if __name__=="__main__":

    r = reader("./mesh.su2")
    triang = mtri.Triangulation(r.x, r.y, r.elem)

    plt.figure()    
    plt.triplot(triang, 'ko-') 
    plt.axis('equal') 
    plt.show()


