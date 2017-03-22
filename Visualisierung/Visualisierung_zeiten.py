import matplotlib.pyplot as pyplot

python = open("Visualisierung_Daten/ergebnis_julia_matrixfrei.txt", "r")
python2 = open("Visualisierung_Daten/ergebnis_julia_matrixfrei0.5.0.txt", "r")
N = []
zeit1 = []
zeit2 = []
for line in python:
    array = line.split(";")
    N.append(array[0])
    zeit1.append(array[4])
i=0
for line2 in python2:
    array2 = line2.split(";")
    zeit2.append((float(zeit1[i])-float(array2[4]))/float(zeit1[i]))
    i=i+1
python.close()
python2.close()

python = open("Visualisierung_Daten/ergebnis_julia_matrixfrei.txt", "r")
python2 = open("Visualisierung_Daten/ergebnis_julia_matrixfrei0.5.0.txt", "r")
N = []
zeit1 = []
zeit2 = []
for line in python:
    array = line.split(";")
    N.append(array[0])
    zeit1.append(array[4])
i=0
for line2 in python2:
    array2 = line2.split(";")
    zeit2.append((float(zeit1[i])-float(array2[4]))/float(zeit1[i]))
    i=i+1
python.close()
python2.close()

pyplot.ylabel("Faktor")
pyplot.xlabel("N")

pyplot.semilogx(N,zeit2,"gs-", linewidth=1, markersize=7, label='(0.4.8-0.5.0)/0.4.8')

pyplot.legend(numpoints=1,loc=1)
pyplot.grid()

pyplot.tight_layout()
pyplot.savefig('Visualisierung_fertig/jacobi_matrixfrei_zeitenfaktor2.png', bbox_inches='tight')