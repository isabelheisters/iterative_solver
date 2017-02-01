import matplotlib.pyplot as pyplot
import matplotlib.patches as patches

#Der Eintrag in die Datei sieht folgendermassen aus:
#k; Dimesionen; Residuum; Verfahren; Iterationen; Zeit; Fehler; Spektralradius

python = open("ergebnisjacobi_python_neu.txt", "r")
xpy11 = []
ypy11 = []
xpy12 = []
ypy12 = []

for line in python:
    array = line.split(";")
    if (float(array[0]) - (float(array[1]) / 4) < 1):
        xpy11.append(array[4])
        ypy11.append(array[2])
    elif (float(array[0]) - (float(array[1]) / 4*3) < 1):
        xpy12.append(array[4])
        ypy12.append(array[2])

pyplot.ylabel("Residuum")
pyplot.xlabel("Iterationen")

pyplot.semilogy(xpy11,ypy11,"bo-", linewidth=2, markersize=10, label='400 k=3/4')
pyplot.semilogy(xpy12, ypy12,"go-", linewidth=2, markersize=10, label='1200 k=1/4')

pyplot.legend(numpoints=1)
pyplot.show()