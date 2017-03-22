import matplotlib.pyplot as pyplot
import numpy as np

python = open("Visualisierung_Daten/ergebnis_mehrgitter_zeit_fourier.txt", "r")
python2 = open("Visualisierung_Daten/ergebnis_julia_2diterativerloeser_verglaich.txt", "r")

n=[]
zeitr=[]
zeiti=[]
n2=[]
zeit=[]

for line in python:
    array = line.split(";")
    print(array)
    ja=[array[3]=="rek", array[3]=="itv"]
    if ja[0]:
        n.append(int(array[4]))
        zeitr.append(float(array[2]))
    elif ja[1]:
        zeiti.append(float(array[2]))

for line in python2:
    array = line.split(";")
    ja=[array[5]=="RGBS", array[5]=="NGS", array[5]=="OJC"]
    if ja[2]:
        n2.append(int(array[1]))
        zeit.append(float(array[2]))

python.close()
python2.close()
pyplot.ylabel("Anzahl Iterationen")
pyplot.xlabel("k")

pyplot.plot(n,zeitr,"gv-", linewidth=1, markersize=7, label='Mehrgitter')
#pyplot.plot(n,zeiti,"bv-", linewidth=1, markersize=7, label='Mehrgitter iterativ')
pyplot.plot(n2,zeit,"yv-", linewidth=1, markersize=7, label='Jacobi')

pyplot.legend(numpoints=1,loc=1)
pyplot.grid()

pyplot.tight_layout()
pyplot.savefig('Visualisierung_fertig/visualisierung_mehrgitter_gegen_jacobi_iterationen.png', bbox_inches='tight')