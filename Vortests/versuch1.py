import math
from time import *

def maximum(h, x, y):
    ableitungx = "-2*x"
    ableitungy = "-2*y"
    while (x>h  or x<-h) and (y>h or y<-h):
        gradx = eval(ableitungx)
        grady = eval(ableitungy)
        betraggrad = math.sqrt((gradx*gradx)+(grady*grady))
        gradx = (1/betraggrad)*gradx
        grady = (1/betraggrad)*grady
        x = x + gradx*h
        y = y + grady*h
    print("Maximum(", x,"|",y,")")

def ableitungx(funktion, h, x, y):
    betrag = eval(funktion)
    x = x + h
    return (eval(funktion)-betrag)/h

def ableitungy(funktion, h, x,y):
    betrag = eval(funktion)
    y = y + h
    return (eval(funktion) - betrag) / h

def maximum_mit_ableitung(funktion, h, x,y):
    while (x > h or x < -h) and (y > h or y < -h):
        gradx = ableitungx(funktion, h, x, y)
        grady = ableitungy(funktion, h, x, y)
        betraggrad = math.sqrt((gradx * gradx) + (grady * grady))
        gradx = (1 / betraggrad) * gradx
        grady = (1 / betraggrad) * grady
        x = x + gradx * h
        y = y + grady * h
    print("Maximum(", x, "|", y, ")")

t1 = clock()
maximum(0.000001, 2, 3)
t2 = clock()
print("Zeit:", t2-t1,"s")
maximum_mit_ableitung("1-x*x-y*y", 0.000001, 2, 3)
t3 = clock()
print("Zeit:", t3-t2,"s")