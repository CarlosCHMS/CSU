
import sys

class reader():

    def __init__(self, fileName):
        
        ff = open(fileName, "r")
        self.elem = []
        self.p = []
        self.markers = []
    
        jj = 0
    
        for row in ff:
            if row[0:5] == "NDIME":
                
                self.Ndime = int(row[6:-1])
                mod = 0
                ii = 0

            elif row[0:5] == "NELEM":
                
                self.Nelem = int(row[6:-1])
                mod = 1    
                ii = 0
                
            elif row[0:5] == "NPOIN":

                self.Npoin = int(row[6:-1]) 
                mod = 2
                ii = 0

            elif row[0:5] == "NMARK":

                self.Nmark = int(row[6:-1])
                mod = 3
                ii = 0

            else:

                if mod == 1:
                    aux = row.split(' ')
                    if int(aux[0]) == 5:                    
                        self.elem.append([int(aux[1]), int(aux[2]), int(aux[3])])
                    elif int(aux[0]) == 9:                    
                        self.elem.append([int(aux[1]), int(aux[2]), int(aux[3]), int(aux[4])])
                        
                    ii += 1

                elif mod == 2:
                    aux = row.split(' ')
                    self.p.append([float(aux[0]), float(aux[1]), float(aux[2])])
                    ii += 1

                elif mod == 3:
                    if jj == 0:
                        m = marker()
                        m.getTag(row)
                        jj += 1
                    elif jj == 1:
                        m.getNelem(row)
                        jj += 1
                    elif jj >= 2 and jj < 1 + m.Nelem:
                        m.elemAppend(row)
                        jj += 1
                    else:
                        m.elemAppend(row)
                        self.markers.append(m)
                        jj = 0


        ff.close()
        
        self.xy()
      

    def xy(self):
    
        self.x = []
        self.y = []

        for pp in self.p:

            self.x.append(pp[0])
            self.y.append(pp[1])

        return None
    
    def writeFile(self):
    
        Ntab = 2 + len(self.markers)
        
        ff = open('./mesh.csv', 'w')

        ff.write("%i, \n" % (Ntab))
        
        for ii in range(0, len(self.markers)):
            ff.write("%i, %i, %i,\n" % (ii, self.markers[ii].Nelem, 2))
            for jj in range(0, self.markers[ii].Nelem):
                ff.write("%i, %i,\n" % (self.markers[ii].elem[jj][0], self.markers[ii].elem[jj][1]))
        
        ff.write("%i, %i, %i,\n" % (Ntab-2, len(self.elem), 3))
        for jj in range(0, len(self.elem)):
            ff.write("%i, %i, %i,\n" % (self.elem[jj][0], self.elem[jj][1], self.elem[jj][2]))

        ff.write("%i, %i, %i,\n" % (Ntab-1, len(self.p), 2))
        for jj in range(0, len(self.p)):
            ff.write("%.10e, %.10e,\n" % (self.p[jj][0], self.p[jj][1]))

        ff.close()
        
        return None
        
        
class marker():
    
    def __init__(self):

        self.tag = ''
        self.Nelem = 0
        self.elem = []

    def getTag(self, row):

        aux = row.split(' ')
        self.tag = str(aux[1])

        return None

    def getNelem(self, row):

        aux = row.split(' ')
        self.Nelem = int(aux[1])

        return None

    def elemAppend(self, row):

        aux = row.split(' ')
        self.elem.append([int(aux[1]), int(aux[2])])
        
    def getXY(self, mesh):
    
        self.x = []
        self.y = []
        
        for e in self.elem:
            self.x.append(mesh.x[e[0]])
            self.y.append(mesh.y[e[0]])
            
        return None


if __name__=="__main__":

    if len(sys.argv) < 1:
        print("Insert the input file")
        exit(0)

    r = reader(sys.argv[1])

