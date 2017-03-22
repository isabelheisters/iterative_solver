import matplotlib.pyplot as pyplot
import numpy as np

python = open("Visualisierung_Daten/ergebnis_mehrgitter_residuum_fourier.txt", "r")
python2 = open("Visualisierung_Daten/ergebnis_julia_2diterativerloeser_residuum.txt", "r")

n=[]
zeitr=[]
n2=[]
zeit=[]

for line in python:
    array = line.split(";")
    n.append(int(array[0]))
    zeitr.append(float(array[1]))

print(n, zeitr)

for line in python2:
    array = line.split(";")
    n2.append(int(array[0]))
    zeit.append(float(array[1]))

pyplot.ylabel("Residuum")
pyplot.xlabel("Iteration")

pyplot.semilogy(n,zeitr,"gv-", linewidth=1, markersize=7, label='Mehrgitter')
pyplot.semilogy(n2[1:100],zeit[1:100],"bv-", linewidth=1, markersize=7, label='Jacobi')

pyplot.legend(numpoints=1,loc=1)
pyplot.grid()

pyplot.tight_layout()
pyplot.savefig('Visualisierung_fertig/visualisierung_mehrgitter_gegen_jacobi_residuum.png', bbox_inches='tight')