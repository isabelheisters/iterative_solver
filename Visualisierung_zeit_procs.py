import matplotlib.pyplot as pyplot
import numpy as np

python = open("../Paralellisierung/Daten/matrixmultiplikation_zweidimensional_1.txt", "r")
python2 = open("../Paralellisierung/Daten/matrixmultiplikation_zweidimensional_4.txt", "r")
python3 = open("../Paralellisierung/Daten/matrixmultiplikation_zweidimensional_9.txt", "r")
python4 = open("../Paralellisierung/Daten/matrixmultiplikation_zweidimensional_16.txt", "r")

n=[]
zeit1=[]
zeit4=[]
zeit9=[]
zeit16=[]

for line in python:
    array = line.split(";")
    n.append(int(array[0]))
    zeit1.append(float(array[3]))
for line in python2:
    array = line.split(";")
    zeit4.append(float(array[3]))
for line in python3:
    array = line.split(";")
    zeit9.append(float(array[3]))
for line in python4:
    array = line.split(";")
    zeit16.append(float(array[3]))

python.close()
python2.close()
python3.close()
python4.close()

pyplot.ylabel("Zeit")
pyplot.xlabel("N")

pyplot.semilogx(n,zeit1,"gv-", linewidth=1, markersize=7, label='1 Prozessor')
pyplot.semilogx(n,zeit4,"gs-", linewidth=1, markersize=7, label='4 Prozessor')
pyplot.semilogx(n,zeit9,"g+-", linewidth=1, markersize=7, label='9 Prozessor')
pyplot.semilogx(n,zeit16,"gp-", linewidth=1, markersize=7, label='16 Prozessor')

pyplot.legend(numpoints=1,loc=1)
pyplot.grid()

pyplot.tight_layout()
pyplot.savefig('Visualisierung_fertig/visualisierung_matrixmultiplikation_zeit_procs.png', bbox_inches='tight')