# Erzeugt eine Datei mit 5000 Zufallszahlen zwischen 0 und 1
import random

datei = open("zufallszahlen.txt", "w")
i=0
while(i<10000000):
    i=i+1
    datei.write(str(random.random())+ "\n")