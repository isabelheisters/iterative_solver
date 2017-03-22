import matplotlib.pyplot as pyplot

python = open("Visualisierung_Daten/ergebnis_interpolation.txt", "r")
N = []
error = []
hquadrat = []
zeit=[]
for line in python:
    array = line.split(";")
    N.append(array[0])
    error.append(array[1])
    hquadrat.append(array[2])
    zeit.append(array[3])

python.close()

pyplot.ylabel("Zeit")
pyplot.xlabel("N")

#pyplot.loglog(N,error,"gs-", linewidth=1, markersize=7, label='Abstand zwischen Abschaetzung und genauem Ergebnis')
#pyplot.loglog(N,hquadrat,"bs-", linewidth=1, markersize=7, label='h^2')
pyplot.loglog(N,zeit,"bs-", linewidth=1, markersize=7, label='Zeit')

pyplot.legend(numpoints=1,loc=1, fontsize='xx-small')
pyplot.grid()

pyplot.tight_layout()
pyplot.savefig('Visualisierung_fertig/interpolation_zeit.png', bbox_inches='tight')