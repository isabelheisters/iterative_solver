import matplotlib.pyplot as pyplot
import matplotlib.patches as patches

#Der Eintrag in die Datei sieht folgendermassen aus:
#k; Dimesionen; Residuum; Verfahren; Iterationen; Zeit; Fehler; Spektralradius

python = open("ergebnisiterationen.txt", "r")
xgs = []
ygs = []
xsor = []
ysor = []
xjc = []
yjc = []
xojc = []
yojc = []



for line in python:
    array=line.split(";")

    if(array[2]==" GS"):
        xgs.append(int(array[0]))
        ygs.append(array[3])
    elif (array[2] == " SOR"):
        xsor.append(int(array[0]))
        ysor.append(array[3])
    elif (array[2] == " JC"):
        xjc.append(int(array[0]))
        yjc.append(array[3])
    elif (array[2] == " OJC"):
        xojc.append(int(array[0]))
        yojc.append(array[3])
python.close()
#f, (eins, zwei,drei,vier) = pyplot.subplots(1,4, sharex=True)

pyplot.ylabel("Anzahl der Iterationen")
pyplot.xlabel("Dimensionen")

pyplot.loglog(xgs,ygs,"bo-", linewidth=2, markersize=10, label="Gauss Seidel")
#eins.set_title("Gauss-Seidel")

pyplot.loglog(xsor,ysor,"gd-", linewidth=2, markersize=10, label="SOR")
#zwei.set_title("SOR")

pyplot.loglog(xjc,yjc,"rs-", linewidth=2, markersize=10, label="Jacobi")
#drei.set_title("Jacobi")

pyplot.loglog(xojc,yojc,"yv-", linewidth=2, markersize=10, label="Gewichteter Jacobi")
#vier.set_title("Gewichteter Jacobi")

pyplot.legend(numpoints=1,loc=4)
pyplot.grid()

pyplot.tight_layout()
pyplot.savefig('iterationen.png', transparent=True, bbox_inches='tight')
# pyplot.show()