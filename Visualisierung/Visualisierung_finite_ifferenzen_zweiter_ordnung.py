import matplotlib.pyplot as pyplot

python = open("Visualisierung_Daten/ergebnis_finite_differenzen_zweite_ordnung.txt", "r")
N = []
error = []
hquadrat = []
for line in python:
    array = line.split(";")
    N.append(array[0])
    error.append(array[1])
    hquadrat.append(array[2])

python.close()

pyplot.ylabel("Error")
pyplot.xlabel("N")

pyplot.loglog(N,error,"gs-", linewidth=1, markersize=7, label='Abstand zwischen Abschaetzung und genauem Ergebnis')
pyplot.loglog(N,hquadrat,"bs-", linewidth=1, markersize=7, label='h^2')

pyplot.legend(numpoints=1,loc=1, fontsize='xx-small')
pyplot.grid()

pyplot.tight_layout()
pyplot.savefig('Visualisierung_fertig/finite_elemente_zweiter_ordnung_error.png', bbox_inches='tight')