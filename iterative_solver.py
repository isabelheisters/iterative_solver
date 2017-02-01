import numpy as np
from scipy import sparse
from scipy.sparse import linalg as la
import math
import time


def do_iteration_neu(A,b,x,P,koeffizient):
    k = 0
    #L = la.splu(P)
    while (np.linalg.norm(b - A.dot(x),np.inf) >= koeffizient):
        v = la.spsolve(P, b-A.dot(x))
        #v = L.solve(b-A.dot(x))
        x = x + v
        k += 1
    return x, k

def make_matrix(dimension):
    vorfaktor = (dimension+1)*(dimension+1)
    stencil = [vorfaktor*-1, vorfaktor*2, vorfaktor*-1]
    A = sparse.diags(stencil, [-1, 0, 1], shape=(dimension, dimension))
    return sparse.csc_matrix(A)


def richardson_P(A):
    return sparse.eye(A.shape[0])

def jacobi_P(A):
    return sparse.diags(A.diagonal(),0,format='csc')

def jacobi_P_mit_w(A):
    w = 2 / 3
    return 1/w * sparse.diags(A.diagonal(),0,format='csc')

def gaussseidel_P(A):
    return sparse.csc_matrix(sparse.tril(A, 0))

def gaussseidel_P_mit_w(A):
    #w = 2/(1+math.sin(math.pi*(1/(A.shape[0]+1))))
    w=1.0#2/(1+math.sin(1.0/(n+1))) #geaendert fuer fouriermode
    return (1/w)*sparse.diags(A.diagonal(),0,format='csc')+sparse.csc_matrix(sparse.tril(A, -1))

def spektralradius(A,P):
    Binv = sparse.linalg.inv(P)
    ergebnis = sparse.eye(A.shape[0])-Binv.dot(A)
    return max(abs(np.linalg.eigvals(ergebnis.todense())))

def durschnittlicher_fehler(A,b,x):
    return np.linalg.norm(x,np.inf)
    # return np.linalg.norm(x-(la.spsolve(A,b)),np.inf)

def mach_was(x, residuum, art):
    print(x.shape[0])
    A = make_matrix(x.shape[0])
    b = np.zeros(x.shape[0])
    if(art == "GS"):
        datei.write("GS; ")
        P = gaussseidel_P(A)
    elif(art == "SOR"):
        datei.write("SOR; ")
        P = gaussseidel_P_mit_w(A)
    elif(art == "JC"):
        datei.write("JC; ")
        P = jacobi_P(A)
    elif(art == "OJC"):
        datei.write("OJC; ")
        P = jacobi_P_mit_w(A)
        #print(eigenwert(A, P,5/10))
        #1/10 hier hart eingecodet wegen k=1/10
    elif(art =="RI"):
        datei.write("RI; ")
        P = richardson_P(A)
    else:
        datei.write("UNGUELTIG; ")
        return 0
    l=0

    zeiten = np.zeros(25)
    k = 0
    xerg =0
    while(l<len(zeiten)):
        t1 = time.time()
        xerg, k = do_iteration_neu(A,b,x,P,residuum)
        t2 = time.time()
        zeiten[l]=t2-t1
        l = l+1
    datei.write(str(k)+"; ")
    durchschnitt = 0
    for element in zeiten:
        durchschnitt = durchschnitt + element
    durchschnitt = durchschnitt/len(zeiten)
    datei.write(str(durchschnitt)+"; ")
    fehler = durschnittlicher_fehler(A, b, xerg)
    u = spektralradius(A, P)
    datei.write(str(fehler)+"; ")
    datei.write(str(u)+"\n")
    return 1

def fouriermode(n,k):
    x = np.zeros(n)
    i=0
    while i<n:
        x[i] = math.sin((k*math.pi*(i+1))/(n+1))
        i = i+1
    return x

def eigenwert(A, P, k):
    Binv = sparse.linalg.inv(P)
    ergebnis = sparse.eye(A.shape[0]) - Binv.dot(A)
    anzahl = int(A.shape[0] * k)
    val, vec = la.eigs(ergebnis, k=anzahl)
    print(ergebnis.shape)
    print(abs(val))
    exit()

    return abs(val[anzahl ])

#print("Bitte geben sie die Dimension der Matrix an, die Groesse des Residuums und die Art der Iteration. Das bedeutet:\nGS ist Gauss Seidel\nRI ist Richardson\nJC ist Jacobi\nOJC ist gewichteter Jacobi\nSOR ist SOR(Optimierter Gauss Seidel) ")
#abbruch = False
#einlesen = open("zufallszahlen.txt", "r")
#dimensionen = 100
#residuum = 0.1
#x = np.zeros(dimensionen)
#i=0
#while(i<dimensionen):
#    x[i] = float(einlesen.readline())
 #   i = i+1
#i=0
#j=0
#einlesen.close()
# print(x)

#datei = open("beispielpython.txt", "w")

#Der Eintrag in die Datei sieht folgendermassen aus:
#Dimension; Residuum; Verfahren; Zeit; Iterationen; Fehler; Spektralradius
#while(j<5):
#while(k<dimensionen):
    #residuum =0.1
    #einlesen = open("zufallszahlen.txt", "r")
    #x1 = np.zeros(dimensionen)
    #i = 0
    #while (i < dimensionen):
     #   x1[i] = float(einlesen.readline())
     #   i = i + 1
    #i = 0
    #einlesen.close()
    #k = dimensionen / 4
    #x2 = fouriermode(dimensionen, k)
    #k = dimensionen / 10
    #print(k)
    #x2 = fouriermode(dimensionen, k)
    #while(i<1):
        #datei.write(str(dimensionen)+ "; ")
        #datei.write(str(residuum)+"; ")
        #mach_was(x1, residuum, "GS")
        #datei.write(str(k)+"; ")
        #datei.write(str(dimensionen)+ "; ")
        #datei.write(str(residuum) + "; ")
        #mach_was(x1, residuum, "SOR")
        #datei.write(str(dimensionen)+ "; ")
        #datei.write(str(residuum) + "; ")
        #mach_was(x1, residuum, "JC")
        #datei.write(str(dimensionen)+ "; ")
        #datei.write(str(residuum) + "; ")
        #mach_was(x2, residuum, "OJC")
       # i = i+1
        #residuum = residuum/10
    #j=j+1
#k=k+(dimensionen/4)
#    j=j+1
    #dimensionen = dimensionen*2
#datei.close()

#while(abbruch == False):
#    print("Abbrechen?(Fuer Abbruch 0 eingeben) ")
#    examp = input()
#    if(examp == "0"):
#        abbruch = True
#    else:
#        print("\nDimensionen: ")
#        dimensionen = int(input())
#        print("Residuum: ")
#        residuum = float(input())
#        print("Art: ")
#        art = input()
#        b = np.zeros(dimensionen)
#        b[dimensionen-1] = 49
#        mach_was(b,residuum, art)

datei = open("Visualisierung/Visualisierung_Daten/ergebnis_python_normal_ohne_lu.txt", "w")
n=10
k=4
while(n<20000):
    print(n)
    datei.write(str(n)+";")
    datei.write(str(k)+ ";")
    A=make_matrix(n)
    #P2=gaussseidel_P_mit_w(A)
    P=jacobi_P_mit_w(A)
    koeffizient=0.1
    b=np.zeros(n)
    zeiten=np.zeros(25)
    anzahl=0
    for i in range(len(zeiten)):
        x = fouriermode(n, k)
        t1=time.time()
        _, anzahl=do_iteration_neu(A, b, x, P, koeffizient)
        t2=time.time()
        zeiten[i]=t2-t1
    datei.write(str(anzahl)+ ";")
    add=0
    for element in zeiten:
        add += element
    datei.write(str(add/zeiten.shape[0])+ ";")
    datei.write(str((add / zeiten.shape[0])/anzahl)+ ";")
    datei.write("OJC\n")

    #datei.write(str(n) + ";")
    #datei.write(str(k) + ";")
    #zeiten2 = np.zeros(25)
    # print(n)
    #for i in range(len(zeiten2)):
     #   x = fouriermode(n, k)
     #   t1 = time.time()
     #   _, anzahl = do_iteration_neu(A, b, x, P2, koeffizient)
      #  t2 = time.time()
      #  zeiten2[i] = t2 - t1
      #  print('Anzahl: %s' %anzahl)
    #datei.write(str(anzahl) + ";")
    #add = 0
    #for element in zeiten2:
    ##    add += element
    #datei.write(str(add / zeiten2.shape[0]) + ";")
    #datei.write(str((add / zeiten2.shape[0]) / anzahl) + ";")
    #datei.write("SOR\n")


    n*=2
    k*=2

datei.close()