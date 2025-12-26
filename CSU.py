# Programed in Python3 but run in python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 15:39:49 2020

@author: carlos
"""


import sys
import subprocess as sp
import os

if __name__=="__main__":

    if len(sys.argv) < 1:
        print("Insert the input file")
        exit(0)
    
    path = sys.argv[1]
                
    print("\nCSU - A CFD code for unstructured meshs:\n")    

    os.system("rm ./executable %ssolution.csv" % (path))
    os.system("gcc ./readTables.c ./utils.c ./mesh.c ./input.c boundary.c ./solver.c ./flux.c ./sa.c ./saCC.c ./implicit.c ./gasprop.c ./limiter.c ./laminar.c ./sst.c ./sstTrans.c ./main.c -o ./executable -lm -fopenmp -Wall -O3")
    os.system("./executable %s" % path)
    os.system("python3 %sanalisys.py %s" % (path, path))
