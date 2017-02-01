import matplotlib.pyplot as pyplot
import matplotlib.patches as patches

#Der Eintrag in die Datei sieht folgendermassen aus:
#Dimesionen; Residuum; Verfahren; Iterationen; Zeit; Fehler; Spektralradius

python = open("ergebnis_python.txt", "r")

xpy1 = []
ypy1 = []
xpy2 = []
ypy2 = []
xpy3 = []
ypy3 = []
xpy4 = []
ypy4 = []
xpy5 = []
ypy5 = []
for line in python:
    array=line.split(";")
    print(array[1])
    if(array[0]=="100"):
        xpy1.append(array[1])
        ypy1.append(array[4])
    elif(array[0]=="200"):
        xpy2.append(array[1])
        ypy2.append(array[4])
    elif (array[0] == "400"):
        xpy3.append(array[1])
        ypy3.append(array[4])
    elif (array[0] == "800"):
        xpy4.append(array[1])
        ypy4.append(array[4])
    elif (array[0] == "1600"):
        xpy5.append(array[1])
        ypy5.append(array[4])