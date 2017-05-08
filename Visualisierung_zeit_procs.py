import matplotlib.pyplot as pyplot
import numpy as np

python = open("/home/zam/heisters/projects/iterative_loeser/Parallelisierung/Daten/jacobi_zweidimensional1.txt", "r")
python2 = open("/home/zam/heisters/projects/iterative_loeser/Parallelisierung/Daten/jacobi_zweidimensional2.txt", "r")
python3 = open("/home/zam/heisters/projects/iterative_loeser/Parallelisierung/Daten/jacobi_zweidimensional4.txt", "r")
python4 = open("/home/zam/heisters/projects/iterative_loeser/Parallelisierung/Daten/jacobi_zweidimensional8.txt", "r")
python5 = open("/home/zam/heisters/projects/iterative_loeser/Parallelisierung/Daten/jacobi_zweidimensional9.txt", "r")

n=[]
zeit1=[]
zeit4=[]
zeit9=[]
zeit16=[]
zeiti=[]

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
for line in python4:
    array = line.split(";")
    zeiti.append(float(array[3]))

python.close()
python2.close()
python3.close()
python4.close()
python5.close()

pyplot.ylabel("Zeit")
pyplot.xlabel("Dimension")

pyplot.loglog(n,zeit1,"gv-", linewidth=1, markersize=7, label='1 Prozessor')
pyplot.loglog(n,zeit4,"gs-", linewidth=1, markersize=7, label='2 Prozessoren')
pyplot.loglog(n,zeit9,"g+-", linewidth=1, markersize=7, label='4 Prozessoren')
pyplot.loglog(n,zeit16,"gp-", linewidth=1, markersize=7, label='8 Prozessoren')
pyplot.loglog(n,zeit16,"g*-", linewidth=1, markersize=7, label='16 Prozessoren')

pyplot.legend(numpoints=1,loc=4)
pyplot.grid()

pyplot.tight_layout()
pyplot.savefig('visualisierung_jacobi_zeit_procs.png', bbox_inches='tight')