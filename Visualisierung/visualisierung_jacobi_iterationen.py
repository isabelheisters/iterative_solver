import matplotlib.pyplot as pyplot
import matplotlib.patches as patches

#Der Eintrag in die Datei sieht folgendermassen aus:
#k; Dimesionen; Residuum; Verfahren; Iterationen; Zeit; Fehler; Spektralradius

python = open("ergebnisjacobi_python_neu.txt", "r")
xpy11 = []
ypy11 = []
xpy12 = []
ypy12 = []
xpy13 = []
ypy13 = []

xpy21 = []
ypy21 = []
xpy22 = []
ypy22 = []
xpy23 = []
ypy23 = []

xpy31 = []
ypy31 = []
xpy32 = []
ypy32 = []
xpy33 = []
ypy33 = []

xpy41 = []
ypy41 = []
xpy42 = []
ypy42 = []
xpy43 = []
ypy43 = []

xpy51 = []
ypy51 = []
xpy52 = []
ypy52 = []
xpy53 = []
ypy53 = []

for line in python:
    array=line.split(";")

    if(array[1]==" 100"):

        if(float(array[0])-(float(array[1])/4)<1):
            xpy11.append(array[4])
            ypy11.append(array[2])
        elif(float(array[0])-(float(array[1])/2)<1):
            xpy12.append(array[4])
            ypy12.append(array[2])
        elif (float(array[0]) - (float(array[1])/4*3)<1):
            xpy13.append(array[4])
            ypy13.append(array[2])
    elif (array[1] == " 200"):
        if (float(array[0]) - (float(array[1]) / 4) < 1):
            xpy21.append(array[4])
            ypy21.append(array[2])
        elif (float(array[0]) - (float(array[1]) / 2) < 1):
            xpy22.append(array[4])
            ypy22.append(array[2])
        elif (float(array[0]) - (float(array[1]) / 4 * 3) < 1):
            xpy23.append(array[4])
            ypy23.append(array[2])
    elif (array[1] == " 400"):
        if (float(array[0]) - (float(array[1]) / 4) < 1):
            xpy31.append(array[4])
            ypy31.append(array[2])
        elif (float(array[0]) - (float(array[1]) / 2) < 1):
            xpy32.append(array[4])
            ypy32.append(array[2])
        elif (float(array[0]) - (float(array[1]) / 4 * 3) < 1):
            xpy33.append(array[4])
            ypy33.append(array[2])
    elif (array[1] == " 800"):
        if (float(array[0]) - (float(array[1]) / 4) < 1):
            xpy41.append(array[4])
            ypy41.append(array[2])
        elif (float(array[0]) - (float(array[1]) / 2) < 1):
            xpy42.append(array[4])
            ypy42.append(array[2])
        elif (float(array[0]) - (float(array[1]) / 4 * 3) < 1):
            xpy43.append(array[4])
            ypy43.append(array[2])
    elif (array[1] == " 1600"):
        if (float(array[0]) - (float(array[1]) / 4) < 1):
            xpy51.append(array[4])
            ypy51.append(array[2])
        elif (float(array[0]) - (float(array[1]) / 2) < 1):
            xpy52.append(array[4])
            ypy52.append(array[2])
        elif (float(array[0]) - (float(array[1]) / 4 * 3) < 1):
            xpy53.append(array[4])
            ypy53.append(array[2])
python.close()

f, (eins, zwei,drei,vier,fuenf) = pyplot.subplots(1,5, sharex=True, sharey=True)

pyplot.ylabel("Residuum")
pyplot.xlabel("Iterationen")

eins.semilogy(xpy11,ypy11,"bo-", linewidth=2, markersize=10, label='k=1/4')
eins.semilogy(xpy12,ypy12,"go-", linewidth=2, markersize=10, label='k=1/2')
eins.semilogy(xpy13,ypy13,"ro-", linewidth=2, markersize=10, label='k=3/4')
eins.set_title("100")

zwei.semilogy(xpy21,ypy21,"bo-", linewidth=2, markersize=10, label='k=1/4')
zwei.semilogy(xpy22,ypy22,"go-", linewidth=2, markersize=10, label='k=1/2')
zwei.semilogy(xpy23,ypy23,"ro-", linewidth=2, markersize=10, label='k=3/4')
zwei.set_title("200")

drei.semilogy(xpy31,ypy31,"bo-", linewidth=2, markersize=10, label='k=1/4')
drei.semilogy(xpy32,ypy32,"go-", linewidth=2, markersize=10, label='k=1/2')
drei.semilogy(xpy33,ypy33,"ro-", linewidth=2, markersize=10, label='k=3/4')
drei.set_title("400")

vier.semilogy(xpy41,ypy41,"bo-", linewidth=2, markersize=10, label='k=1/4')
vier.semilogy(xpy42,ypy42,"go-", linewidth=2, markersize=10, label='k=1/2')
vier.semilogy(xpy43,ypy43,"ro-", linewidth=2, markersize=10, label='k=3/4')
vier.set_title("800")

fuenf.semilogy(xpy51,ypy51,"bo-", linewidth=2, markersize=10, label='k=1/4')
fuenf.semilogy(xpy52,ypy52,"go-", linewidth=2, markersize=10, label='k=1/2')
fuenf.semilogy(xpy53,ypy53,"ro-", linewidth=2, markersize=10, label='k=3/4')
fuenf.set_title("1600")

pyplot.legend(numpoints=1)
pyplot.show()