import matplotlib.pyplot as pyplot
import numpy as np

python = open("Visualisierung_Daten/ergebnis_mehrgitter_zeit.txt", "r")
python2 = open("Visualisierung_Daten/ergebnis_mehrgitter_objektorientiert_zeit.txt", "r")
python3 = open("Visualisierung_Daten/ergebnis_mehrgitter_zeit_rbgs.txt", "r")
python4 = open("Visualisierung_Daten/ergebnis_mehrgitter_objektorientiert_zeit_rbgs.txt", "r")

n=[]
zeitr=[]
zeiti=[]
zeitroo=[]
zeitioo=[]
zeitrgs=[]
zeitigs=[]
zeitroogs=[]
zeitioogs=[]
for line in python:
    array = line.split(";")
    ja=[array[3]=="rek", array[3]=="itv"]
    if ja[0]:
        n.append(int(array[0]))
        zeitr.append(float(array[1]))
    elif ja[1]:
        zeiti.append(float(array[1]))
for line in python2:
    array = line.split(";")
    print(array)
    ja=[array[3]=="rek", array[3]=="itv"]
    if ja[0]:
        zeitroo.append(float(array[1]))
    elif ja[1]:
        zeitioo.append(float(array[1]))
for line in python3:
    array = line.split(";")
    ja=[array[3]=="rek", array[3]=="itv"]
    if ja[0]:
        zeitrgs.append(float(array[1]))
    elif ja[1]:
        zeitigs.append(float(array[1]))
for line in python4:
    array = line.split(";")
    print(array)
    ja=[array[3]=="rek", array[3]=="itv"]
    if ja[0]:
        zeitroogs.append(float(array[1]))
    elif ja[1]:
        zeitioogs.append(float(array[1]))

zeitroo = np.divide(zeitroo,zeiti)
zeitioo = np.divide(zeitioo,zeiti)
zeitr = np.divide(zeitr,zeiti)
zeitroogs = np.divide(zeitroogs,zeiti)
zeitioogs = np.divide(zeitioogs,zeiti)
zeitrgs = np.divide(zeitrgs,zeiti)
zeitigs = np.divide(zeitigs,zeiti)
end=len(n)
anf=1
python.close()
pyplot.ylabel("Prozentual(/iterativ)")
pyplot.xlabel("N")

pyplot.semilogx(n[anf:end],zeitroo[anf:end],"gv-", linewidth=1, markersize=7, label='rekursiv OO')
pyplot.semilogx(n[anf:end],zeitioo[anf:end],"bv-", linewidth=1, markersize=7, label='iterativ OO')
pyplot.semilogx(n[anf:end],zeitr[anf:end],"gs-", linewidth=1, markersize=7, label='rekursiv')
pyplot.semilogx(n[anf:end],zeitroogs[anf:end],"g+-", linewidth=1, markersize=7, label='rekursiv OO Red Black Gauss Seidel')
pyplot.semilogx(n[anf:end],zeitioogs[anf:end],"b+-", linewidth=1, markersize=7, label='iterativ OO Red Black Gauss Seidel')
#pyplot.semilogx(n[anf:end],zeitrgs[anf:end],"gp-", linewidth=1, markersize=7, label='rekursiv Red Black Gauss Seidel')
pyplot.semilogx(n[anf:end],zeitigs[anf:end],"bp-", linewidth=1, markersize=7, label='iterativ Red Black Gauss Seidel')
#pyplot.loglog(n[anf:end],zeiti[anf:end],"bs-", linewidth=1, markersize=7, label='iterativ')


pyplot.legend(numpoints=1,loc=1)
pyplot.grid()

pyplot.tight_layout()
pyplot.savefig('Visualisierung_fertig/visualisierung_mehrgitter_zeit_prozentual2.png', bbox_inches='tight')