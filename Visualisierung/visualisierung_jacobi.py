import matplotlib.pyplot as pyplot
import matplotlib.patches as patches

#Der Eintrag in die Datei sieht folgendermassen aus:
#Dimesionen; Residuum; Verfahren; Iterationen; Zeit; Fehler; Spektralradius

python = open("ergebnisjacobi_pythonzeit.txt", "r")
xpy1 = []
ypy1 = []
xpy2 = []
ypy2 = []
xpy3 = []
ypy3 = []
xpy4 = []
ypy4 = []
xpy5 = []
ypy5 = []
for line in python:
    array=line.split(";")
    print(array[1])
    if(array[1]==" 100"):
        xpy1.append(array[2])
        ypy1.append(array[5])
    elif(array[1]==" 200"):
        xpy2.append(array[2])
        ypy2.append(array[5])
    elif (array[1] == " 400"):
        xpy3.append(array[2])
        ypy3.append(array[5])
    elif (array[1] == " 800"):
        xpy4.append(array[2])
        ypy4.append(array[5])
    elif (array[1] == " 1600"):
        xpy5.append(array[2])
        ypy5.append(array[5])

python = open("ergebnisjacobi_pythonzeit_neu.txt", "r")
xpy1neu = []
ypy1neu = []
xpy2neu = []
ypy2neu = []
xpy3neu = []
ypy3neu = []
xpy4neu = []
ypy4neu = []
xpy5neu = []
ypy5neu = []
for line in python:
    array=line.split(";")
    if(array[1]==" 100"):
        xpy1neu.append(array[2])
        ypy1neu.append(array[5])
    elif(array[1]==" 200"):
        xpy2neu.append(array[2])
        ypy2neu.append(array[5])
    elif (array[1] == " 400"):
        xpy3neu.append(array[2])
        ypy3neu.append(array[5])
    elif (array[1] == " 800"):
        xpy4neu.append(array[2])
        ypy4neu.append(array[5])
    elif (array[1] == " 1600"):
        xpy5neu.append(array[2])
        ypy5neu.append(array[5])

python.close()
julia = open("ergebnisjacobi_juliazeit.txt", "r")
xjl1 = []
yjl1= []
xjl2 = []
yjl2= []
xjl3 = []
yjl3= []
xjl4 = []
yjl4= []
xjl5 = []
yjl5= []
for line in julia:
    print(array)
    array=line.split(";")
    if (array[1] == " 100"):
        xjl1.append(array[2])
        yjl1.append(array[5])
    elif (array[1] == " 200"):
        xjl2.append(array[2])
        yjl2.append(array[5])
    elif (array[1] == " 400"):
        xjl3.append(array[2])
        yjl3.append(array[5])
    elif (array[1] == " 800"):
        xjl4.append(array[2])
        yjl4.append(array[5])
    elif (array[1] == " 1600"):
        xjl5.append(array[2])
        yjl5.append(array[5])
julia.close()

julia = open("ergebnisjacobi_juliazeit_neu.txt", "r")
xjlneu1 = []
yjlneu1= []
xjlneu2 = []
yjlneu2= []
xjlneu3 = []
yjlneu3= []
xjlneu4 = []
yjlneu4= []
xjlneu5 = []
yjlneu5= []
for line in julia:
    array=line.split(";")
    if (array[1] == " 100"):
        xjlneu1.append(array[2])
        yjlneu1.append(array[5])
    elif (array[1] == " 200"):
        xjlneu2.append(array[2])
        yjlneu2.append(array[5])
    elif (array[1] == " 400"):
        xjlneu3.append(array[2])
        yjlneu3.append(array[5])
    elif (array[1] == " 800"):
        xjlneu4.append(array[2])
        yjlneu4.append(array[5])
    elif (array[1] == " 1600"):
        xjlneu5.append(array[2])
        yjlneu5.append(array[5])
julia.close()

pyplot.ylabel("Zeit")
pyplot.xlabel("Residuum")
#pyplot.semilogx(xpy1,ypy1,"bo-", linewidth=2, markersize=10, label='Python 100')
#pyplot.semilogx(xjl1, yjl1,"go-", linewidth=2, markersize=10, label='Julia 100')
#pyplot.semilogx(xpy1neu,ypy1neu,"ro-", linewidth=2, markersize=10, label='Python 100 neu')
#pyplot.semilogx(xjlneu1,yjlneu1,"yo-", linewidth=2, markersize=10, label='Julia 100 neu')

pyplot.semilogx(xpy2,ypy2,"bs-", linewidth=1, markersize=7, label='Python 200')
pyplot.semilogx(xjl2, yjl2,"gs-", linewidth=1, markersize=7, label='Julia 200')
pyplot.semilogx(xpy2neu,ypy2neu,"rs-", linewidth=1, markersize=7, label='Python 200 neu')
pyplot.semilogx(xjlneu2, yjlneu2,"ys-", linewidth=1, markersize=7, label='Julia 200 neu')

pyplot.semilogx(xpy3,ypy3,"bv-", linewidth=1, markersize=7, label='Python 400')
pyplot.semilogx(xjl3, yjl3,"gv-", linewidth=1, markersize=7, label='Julia 400')
pyplot.semilogx(xpy3neu,ypy3neu,"rv-", linewidth=1, markersize=7, label='Python 400 neu')
pyplot.semilogx(xjlneu3, yjlneu3,"yv-", linewidth=1, markersize=7, label='Julia 400 neu')

pyplot.semilogx(xpy4,ypy4,"b*-", linewidth=1, markersize=7, label='Python 800')
pyplot.semilogx(xjl4, yjl4,"g*-", linewidth=1, markersize=7, label='Julia 800')
pyplot.semilogx(xpy4neu,ypy4neu,"r*-", linewidth=1, markersize=7, label='Python 800 neu')
pyplot.semilogx(xjlneu4, yjlneu4,"y*-", linewidth=1, markersize=7, label='Julia 800 neu')

#pyplot.semilogx(xpy5,ypy5,"bx-", linewidth=1, markersize=7, label='Python 1600')
#pyplot.semilogx(xjl5, yjl5,"gx-", linewidth=1, markersize=7, label='Julia 1600')
pyplot.semilogx(xpy5neu,ypy5neu,"rx-", linewidth=1, markersize=7, label='Python 1600 neu')
pyplot.semilogx(xjlneu5, yjlneu5,"yx-", linewidth=1, markersize=7, label='Julia 1600 neu')

# green_patch = patches.Patch(color='green', label='Julia')
# blue_patch = patches.Patch(color='blue', label='Python')
#hundret_patch = patches.Patch(marker='o', label='100')
#twohundret_patch = patches.Patch(marker='s', label='200')
#threehundret_patch = patches.Patch(marker='v', label='400')
pyplot.legend(numpoints=1,loc=1, fontsize='xx-small')
pyplot.grid()

pyplot.tight_layout()
pyplot.savefig('Jacobineu.png', transparent=True, bbox_inches='tight')