import matplotlib.pyplot as pyplot
import matplotlib.patches as patches

#Der Eintrag in die Datei sieht folgendermassen aus:
#Dimesionen; Residuum; Verfahren; Iterationen; Zeit; Fehler; Spektralradius

python = open("ergebnis_python_normal.txt", "r")
xpyjcn = []
ypyjcn = []
xpysorn = []
ypysorn = []
for line in python:
    array = line.split(";")
    if(array[5]=="OJC"):
        xpyjcn.append(array[0])
        ypyjcn.append(array[4])
    elif(array[5] == "SOR"):
        xpysorn.append(array[0])
        ypysorn.append(array[4])

python = open("ergebnis_python_matrixfrei.txt", "r")
xpyjcm = []
ypyjcm = []
xpysorm = []
ypysorm = []
for line in python:
    array = line.split(";")
    if(array[5]=="OJC"):
        xpyjcm.append(array[0])
        ypyjcm.append(array[4])
    elif(array[5] == "SOR"):
        xpysorm.append(array[0])
        ypysorm.append(array[4])

python = open("ergebnis_julia_normal.txt", "r")
xjljcn = []
yjljcn = []
xjlsorn = []
yjlsorn = []
for line in python:
    array = line.split(";")
    if(array[5]=="OJC"):
        xjljcn.append(array[0])
        yjljcn.append(array[4])
    elif(array[5] == "SOR"):
        xjlsorn.append(array[0])
        yjlsorn.append(array[4])

python = open("ergebnis_julia_matrixfrei.txt", "r")
xjljcm = []
yjljcm = []
xjlsorm = []
yjlsorm = []
for line in python:
    array = line.split(";")
    if(array[5]=="OJC"):
        xjljcm.append(array[0])
        yjljcm.append(array[4])
    elif(array[5] == "SOR"):
        xjlsorm.append(array[0])
        yjlsorm.append(array[4])



pyplot.ylabel("Zeit pro Iteration")
pyplot.xlabel("Dimensionen")

pyplot.semilogx(xpyjcn,ypyjcn,"bs-", linewidth=1, markersize=7, label='Python Jacobi normal')
pyplot.semilogx(xpysorn,ypysorn,"bv-", linewidth=1, markersize=7, label='Python SOR normal')
pyplot.semilogx(xpyjcm,ypyjcm,"b*-", linewidth=1, markersize=7, label='Python Jacobi matrixfrei')
pyplot.semilogx(xpysorm,ypysorm,"bo-", linewidth=1, markersize=7, label='Python SOR matrixfrei')

pyplot.semilogx(xjljcn,yjljcn,"gs-", linewidth=1, markersize=7, label='Julia Jacobi normal')
pyplot.semilogx(xjlsorn,yjlsorn,"gv-", linewidth=1, markersize=7, label='Julia SOR normal')
pyplot.semilogx(xjljcm,yjljcm,"g*-", linewidth=1, markersize=7, label='Julia Jacobi matrixfrei')
pyplot.semilogx(xjlsorm,yjlsorm,"go-", linewidth=1, markersize=7, label='Julia SOR matrixfrei')


pyplot.legend(numpoints=1,loc=1, fontsize='xx-small')
pyplot.grid()

pyplot.tight_layout()
pyplot.savefig('visualisierung.png', transparent=True, bbox_inches='tight')