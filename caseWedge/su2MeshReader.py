


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
                    self.elem.append([int(aux[1]), int(aux[2]), int(aux[3])])
                    ii += 1

                elif mod == 2:
                    aux = row.split(' ')
                    self.p.append([float(aux[0]), float(aux[1]), float(aux[2])])
                    ii += 1

                elif mod == 3:
                    if jj == 0:
                        print('1', row)
                        m = marker()
                        m.getTag(row)
                        jj += 1
                    elif jj == 1:
                        print('2', row)
                        m.getNelem(row)
                        jj += 1
                    elif jj >= 2 and jj < 1 + m.Nelem:
                        print('3', row)
                        m.elemAppend(row)
                        jj += 1
                    else:
                        print('4', row)
                        m.elemAppend(row)
                        self.markers.append(m)
                        jj = 0


                        
                        
                
                        


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


if __name__=="__main__":

    r = reader("./testesu2")

    print(r.Ndime)
    print(r.Nelem)
    print(r.elem)
    print(r.Npoin)
    print(r.p)
    print(r.Nmark)
    for m in r.markers:
        print(m.elem)

