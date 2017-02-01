import copy
import numpy as np
import math
import time

def bedingung(A,b,x):
    n = x.shape[0]
    ergebnis=np.zeros(n)
    ergebnis[0] = A[1] * x[0] + A[2] * x[1]
    for i in range(1,n-1):
        ergebnis[i]=A[0]*x[i-1]+A[1]*x[i]+A[2]*x[i+1]
    ergebnis[-1] = A[0] * x[-2] + A[1] * x[-1]
    return ergebnis

def iterativer_loeser(a,b,x,residuum,w):
    n = x.shape[0]
    anzahl=0
    while(residuum<=np.linalg.norm(b-bedingung(a,b,x),np.inf)):
        xalt=copy.deepcopy(x)
        x[0] = (1 - w) * xalt[0] + ((w / a[1]) * (b[0] - (a[2] * xalt[1])))
        for i in range(1, n - 1):
            x[i]= (1-w)*xalt[i]+((w/a[1])*(b[i]-(a[0]*xalt[i-1]+a[2]*xalt[i+1])))
        x[-1] = (1 - w) * xalt[-1] + ((w / a[1]) * (b[-1] - (a[0] * xalt[-2])))
        anzahl=anzahl+1
    return x,anzahl

def iterativer_loeser_sor(a,b,x,residuum,w):
    n = x.shape[0]
    anzahl=0
    while(residuum<=np.linalg.norm(b-bedingung(a,b,x),np.inf)):
        x[0] = (1 - w) * x[0] + ((w / a[1]) * (b[0] - (a[2] * x[1])))
        for i in range(1, n - 1):
            x[i]= (1-w)*x[i]+((w/a[1])*(b[i]-(a[0]*x[i-1]+a[2]*x[i+1])))
        x[-1] = (1 - w) * x[-1] + ((w / a[1]) * (b[-1] - (a[0] * x[-2])))
        anzahl=anzahl+1
    return x,anzahl

def make_A(n):
    return [-((n + 1) * (n + 1)), 2 * ((n + 1) * (n + 1)), -((n + 1) * (n + 1))]

def fouriermode(n,k):
    x = np.zeros(n)
    i=0
    while i<n:
        x[i] = math.sin((k*math.pi*(i+1))/(n+1))
        i = i+1
    return x

def jacobi(n,k, koeffizient):
   w=2/3
   A = make_A(n)
   b = np.zeros(n)
   zeiten = np.zeros(25)
   for i in range(len(zeiten)):
        x = fouriermode(n, k)
        t1=time.time()
        _, anzahl = iterativer_loeser(A, b,x,koeffizient,w)
        t2=time.time()
        zeiten[i]=(t2-t1)
   datei.write(str(anzahl) + ";")
   add = 0
   for element in zeiten:
      add += element
   datei.write(str(add / zeiten.shape[0]) + ";")
   datei.write(str((add / zeiten.shape[0]) / anzahl) + ";")
   datei.write("OJC\n")

def sor(n,k, koeffizient):
    w=1
    A = make_A(n)
    b = np.zeros(n)
    zeiten = np.zeros(25)
    for i in range(len(zeiten)):
        x = fouriermode(n, k)
        t1 = time.time()
        _, anzahl = iterativer_loeser_sor(A, b, x, koeffizient, w)
        t2 = time.time()
        zeiten[i] = (t2 - t1)
    print('Anzahl: %s' % anzahl)
    datei.write(str(anzahl) + ";")
    add = 0
    for element in zeiten:
        add += element
    datei.write(str(add / zeiten.shape[0]) + ";")
    datei.write(str((add / zeiten.shape[0]) / anzahl) + ";")
    datei.write("SOR\n")


datei = open("Visualisierung/Visualisierung_Daten/ergebnis_python_matrixfrei.txt", "w")
n=10
k=4
koeffizient=0.1
while(n<20000):
    print(n)
    datei.write(str(n)+";")
    datei.write(str(k)+ ";")

    #Hier muss entschieden werden, ob man SOR oder Jacobi oder beides ausfuehren moechte
    #sor(n,k,koeffizient)
    jacobi(n, k, koeffizient)

    n *= 2
    k *= 2

datei.close()