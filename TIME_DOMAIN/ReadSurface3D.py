# ****************************************************************************
#
# Function for surface reading in the advanced time approach of the 
# FW-H method.
# Author: Jaime Vaquero
# September 2017
#
# ****************************************************************************
import numpy as np

##################################################################
def ReadSurface3D(x,y,z,variables,Cartesian=True):
    Ni = len(x)/3
    if (type(variables) is not list):
        variablesList = []
        varibalesList.append(variables)
        variables = variablesList
    variables.append(x)
    variables.append(y)
    variables.append(z)

    Ncol = len(variables)
    if (Cartesian == True):
        print(' ReadSurface3D: using cartesian coordinates')
        outputSortedZ  = SortData(variables,z)
        outputSortedZXY = []
        ilast = -1
        while (ilast<Ni*3-1):
            print(ilast,Ni)
            tempSortedZ,ilast    = TemporalListOfArrays(outputSortedZ,ilast+1)
            print(ilast,Ni)
            tempSortedZX   = SortData(tempSortedZ ,tempSortedZ [-2])
            tempSortedZXY  = SortData(tempSortedZX,tempSortedZX[-3])
            for s in tempSortedZXY:
                outputSortedZXY.append(s)
        finalList = []
        for j in range(0,Ncol):
            finalList.append(outputSortedZXY[j::Ncol])
        print(finalList)    
    else: 
        print(' ReadSurface3D: using cylindrical coordinates')
        outputSortedZ  = SortData(variables,z)

    return outputSortedZXY


def SortData(variables,sortingVariable):
    sortindex  = np.argsort(sortingVariable)
    outputList = [] 
    for x in variables:
       x2 = np.zeros(len(x))
       i  = -1
       for m in sortindex:
           i     = i + 1
           x2[i] = x[m] 
       outputList.append(x2) 
    return outputList

def TemporalListOfArrays(OriginalList,startindex):
    i = startindex
    while (i < len(OriginalList[-1]) and OriginalList[-1][i] == OriginalList[-1][startindex]):
        i = i + 1
    endindex = i-1
    TemporalList = []
    for x in OriginalList:
        x2 = np.zeros(endindex+1-startindex)
        for i in range(startindex,endindex+1):
            x2[i-startindex] = x[i]
        TemporalList.append(x2)

    return TemporalList,endindex
############################################################

file1 = open('Libro1.txt','r')
file1r = file1.read().split()
file1.close()
data = np.zeros(len(file1r))

for i in range(0,len(file1r)):
    data[i] = float(file1r[i])

x = np.zeros(27)
y = np.zeros(27)
z = np.zeros(27)
u = np.zeros(27)
v = np.zeros(27)

x = data[::5 ]
y = data[1::5]
z = data[2::5]
u = data[3::5]
v = data[4::5]

variables = []
variables.append(u)
variables.append(v)
jaja = ReadSurface3D(x,y,z,variables)
print('')
print(' Final List: ')
for i in range(len(jaja)):
    print((jaja[i]))
