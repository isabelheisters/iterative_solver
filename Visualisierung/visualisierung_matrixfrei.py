import matplotlib.pyplot as pyplot
import matplotlib.patches as patches

#Der Eintrag in die Datei sieht folgendermassen aus:
#Dimensionen;prozent;iterationen;ZeitGesamt;ZeitproIteration;Verfahren
x=0
y=4

python = open("Visualisierung_Daten/ergebnis_python_normal.txt", "r")
xpyjcn = []
ypyjcn = []
xpysorn = []
ypysorn = []
for line in python:
    array = line.split(";")
    if(array[5]=="OJC\n"):
        xpyjcn.append(array[x])
        ypyjcn.append(array[y])
    elif(array[5] == "SOR\n"):
        xpysorn.append(array[x])
        ypysorn.append(array[y])

python = open("Visualisierung_Daten/ergebnis_python_matrixfrei.txt", "r")
xpyjcm = []
ypyjcm = []
xpysorm = []
ypysorm = []
for line in python:
    array = line.split(";")
    if(array[5]=="OJC\n"):
        xpyjcm.append(array[x])
        ypyjcm.append(array[y])
    elif(array[5] == "SOR\n"):
        xpysorm.append(array[x])
        ypysorm.append(array[y])

python = open("Visualisierung_Daten/ergebnis_python_normal_ohne_lu.txt", "r")
xpyjclu = []
ypyjclu = []
xpysorlu = []
ypysorlu = []
for line in python:
    array = line.split(";")
    if(array[5]=="OJC\n"):
        xpyjclu.append(array[x])
        ypyjclu.append(array[y])
    elif(array[5] == "SOR\n"):
        xpysorlu.append(array[x])
        ypysorlu.append(array[y])


julia = open("Visualisierung_Daten/ergebnis_julia_normal.txt", "r")
xjljcn = []
yjljcn = []
xjlsorn = []
yjlsorn = []
for line in julia:
    array = line.split(";")
    if(array[5]=="OJC\n"):
        xjljcn.append(array[x])
        yjljcn.append(array[y])
    elif(array[5] == "SOR\n"):
        xjlsorn.append(array[x])
        yjlsorn.append(array[y])

julia = open("Visualisierung_Daten/ergebnis_julia_matrixfrei.txt", "r")
xjljcm = []
yjljcm = []
xjlsorm = []
yjlsorm = []
for line in julia:
    array = line.split(";")
    if(array[5]=="OJC\n"):
        xjljcm.append(array[x])
        yjljcm.append(array[y])
    elif(array[5] == "SOR\n"):
        xjlsorm.append(array[x])
        yjlsorm.append(array[y])

julia = open("Visualisierung_Daten/ergebnis_julia_ohne_lu.txt", "r")
xjljclu = []
yjljclu = []
xjlsorlu = []
yjlsorlu = []
for line in julia:
    array = line.split(";")
    if (array[5] == "OJC\n"):
        xjljclu.append(array[x])
        yjljclu.append(array[y])
    elif (array[5] == "SOR\n"):
        xjlsorlu.append(array[x])
        yjlsorlu.append(array[y])



pyplot.ylabel("Zeit pro Iteration")
pyplot.xlabel("Dimensionen")

pyplot.loglog(xpyjcn,ypyjcn,"bs-", linewidth=1, markersize=7, label='Python Jacobi LU')
#pyplot.loglog(xpysorn,ypysorn,"bv-", linewidth=1, markersize=7, label='Python SOR LU')
pyplot.loglog(xpyjclu,ypyjclu,"bp-", linewidth=1, markersize=7, label='Python Jacobi normal')
#pyplot.loglog(xpysorlu,ypysorlu,"b+-", linewidth=1, markersize=7, label='Python SOR normal')
pyplot.loglog(xpyjcm,ypyjcm,"b*-", linewidth=1, markersize=7, label='Python Jacobi matrixfrei')
#pyplot.loglog(xpysorm,ypysorm,"bo-", linewidth=1, markersize=7, label='Python SOR matrixfrei')

pyplot.loglog(xjljcn,yjljcn,"gs-", linewidth=1, markersize=7, label='Julia Jacobi LU')
#pyplot.loglog(xjlsorn,yjlsorn,"gv-", linewidth=1, markersize=7, label='Julia SOR LU')
pyplot.loglog(xjljclu,yjljclu,"gp-", linewidth=1, markersize=7, label='Julia Jacobi normal')
#pyplot.loglog(xjlsorlu,yjlsorlu,"g+-", linewidth=1, markersize=7, label='Julia SOR normal')
pyplot.loglog(xjljcm,yjljcm,"g*-", linewidth=1, markersize=7, label='Julia Jacobi matrixfrei')
#pyplot.loglog(xjlsorm,yjlsorm,"go-", linewidth=1, markersize=7, label='Julia SOR matrixfrei')


pyplot.legend(numpoints=1,loc=2, fontsize='xx-small')
pyplot.grid()

pyplot.tight_layout()
pyplot.savefig('Visualisierung_fertig/visualisierung_jacobi_iteration_per_time_0.4_all_2.png', transparent=False, bbox_inches='tight')