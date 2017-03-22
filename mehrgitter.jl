import Base.+
import Base.-
import Base.abs

#Datentyp um Matrizen abzuspeichern
type Hmatrix
    element::Array{Float64,2}
end

#Datentyp um Stencils abzuspeichern
type Stencil
    element::Array{Int64,1}
end


function +(a::Hmatrix, b::Hmatrix)
    return Hmatrix(a.element+b.element)
end
function -(a::Hmatrix, b::Hmatrix)
    return Hmatrix(a.element-b.element)
end
function abs(a::Hmatrix)
    return(abs(a.element))
end

# Matrixmultiplikation mit einem Stencil im zweidimensionalen Bereich
function matrixmultiplikation(A::Array{Int64,2},un::Array{Float64,2})
    u=fuege_Rand_hinzu(un)
    gr=size(u)[1]
    matrix=zeros(gr-2,gr-2)
    for i in 2:gr-1
        for j in 2:gr-1
            matrix[i-1,j-1]=A[1]*u[i-1,j]+A[2]*u[i,j-1]+A[3]*u[i,j]+A[4]*u[i,j+1]+A[5]*u[i+1,j]
        end
    end
    return matrix
end

#Ich arbeite mit einer geraenderten Matrix bei diesem Verfahren,
#dies macht die Implementierung leichter.
function fuege_Rand_hinzu(u::Array{Float64,2})
    N=size(u)[1]
    matrix=zeros(N+2,N+2)
    matrix[2:N+1,2:N+1] = u
    return matrix
end

#2-Dimensionaler Jacobi Algorithmus
#Schau in die Doku wenn du nicht weisst wies funktioniert!
function jacobi(x::Array{Float64,2},a::Array{Int64,1},b::Array{Float64,2},anzahl::Int64)
    w=2/3
    x=fuege_Rand_hinzu(x)
    b=fuege_Rand_hinzu(b)
    n=size(x)[1]
    z=0
    while(z<anzahl)
        xalt=copy(x)
        for i in 2:n-1
            for j in 2:n-1
                x[i,j]= (1-w)*xalt[i,j]+((w/a[3])*(b[i,j]-(a[1]*xalt[i-1,j]+a[2]*xalt[i,j-1]+a[4]*xalt[i,j+1]+a[5]*xalt[i+1,j])))
            end
        end
        z=z+1
    end
    x=x[2:end-1, 2:end-1]
    return x
end

function redblackgs(x::Array{Float64,2},a::Array{Int64,1},b::Array{Float64,2},anzahl::Int64)
    x=fuege_Rand_hinzu(x)
    b=fuege_Rand_hinzu(b)
    n=size(x)[1]
    z=0
    while(z<anzahl)
        #Erste Schleife
        for i in 2:(n-1)
            for j in collect(2+(i%2):2:n-1)
                x[i,j]= ((1/a[3])*(b[i,j]-(a[1]*x[i-1,j]+a[2]*x[i,j-1]+a[4]*x[i,j+1]+a[5]*x[i+1,j])))
            end
        end
        #Zweite Schleife
        for i in 2:(n-1)
            for j in collect(3-(i%2):2:n-1)
                x[i,j]= ((1/a[3])*(b[i,j]-(a[1]*x[i-1,j]+a[2]*x[i,j-1]+a[4]*x[i,j+1]+a[5]*x[i+1,j])))
            end
        end
        z=z+1
    end
    x=x[2:end-1, 2:end-1]
    return x
end

#Auch hier arbeite ich mit einer geraenderten Matrix
#Die Matrix wird normal restringiert
function restriktion(x::Array{Float64,2})
    groesse = size(x)[1]
    matrix=zeros(convert(Int,(groesse-1)/2),convert(Int,(groesse-1)/2))
    x=fuege_Rand_hinzu(x)
    k=3
    for i in 1:convert(Int,(groesse-1)/2)
        l=3
        for j in 1:convert(Int,(groesse-1)/2)
            matrix[i,j]=0.5*x[k,l]+0.125*x[k-1,l]+0.125*x[k+1,l]+0.125*x[k,l-1]+0.125*x[k,l+1]
            l=l+2
        end
        k=k+2
    end
    return matrix
end

#Die feinere Matrix wird in die groebere Matrix gesetzt
function einsetzen(matrix::Array{Float64,2})
    groesse= size(matrix)[1]
    neue=zeros(groesse*2+3, groesse*2+3)
    k=1
    for i in range(3,2,groesse)
        l=1
        for j in range(3,2,groesse)
            neue[i,j]=matrix[k,l]
            l=l+1
        end
        k=k+1
     end
    return(neue)
end

#Die Matrix wird beim einsetzen geraendert und es erfolgt eine Prolongation
function prolongation(xi::Array{Float64,2})
    x=einsetzen(xi)
    groesse= size(x)[1]
    #1. Zeilen in denen schon Werte stehen werden bearbeitet
    for i in range(3,2,convert(Int,(groesse-3)/2))
        for j in range(2,2,convert(Int,(groesse-3)/2+1))
            x[i,j]=0.5*x[i, j-1]+0.5*x[i, j+1]
        end
    end
    #2. Zeilen in denen noch keine Werte stehen werden bearbeitet
    for i in range(2,2,convert(Int,(groesse-3)/2+1))
        for j in range(2,1,groesse-2)
            x[i,j]=0.5*x[i-1, j]+0.5*x[i+1, j]
        end
    end
    xi=x[2:end-1, 2:end-1]
    return xi
end

function mehrgitter(A::Array{Int64,2},varray::Array{Array{Float64,2},1}, farray::Array{Array{Float64,2},1},iterationen::Int64, anzahl::Int64)
    #Rekursiver Anker
    if(anzahl==1)
        varray[anzahl][1,1]=farray[anzahl][1,1]/A[anzahl,3]
        return 0
    end
    anzahl=anzahl-1
    
    #Vorglaetten
    v=jacobi(varray[anzahl+1],vec(A[anzahl+1:anzahl+1, 1:end]),farray[anzahl+1],iterationen)
    varray[anzahl+1]=v
    
    #residuum
    r=farray[anzahl+1]-matrixmultiplikation(A[anzahl+1:anzahl+1, 1:end],varray[anzahl+1])
    
    #restriktion
    f= restriktion(r)
    farray[anzahl]=f
    
    #rekursion
    mehrgitter(A,varray, farray,iterationen,anzahl)
    
    #prolongation
    varray[anzahl]=prolongation(varray[anzahl])
    varray[anzahl+1]+=varray[anzahl]
    #Nachglaetten
    varray[anzahl+1]=jacobi(varray[anzahl+1],vec(A[anzahl+1:anzahl+1, 1:end]),farray[anzahl+1],iterationen)
    
end

function mehrgitter_iterativ(A::Array{Int64,2},varray::Array{Array{Float64,2},1}, farray::Array{Array{Float64,2},1},iterationen::Int64, anz::Int64)
    anzahl=anz
    
    while(anzahl>1)
        anzahl=anzahl-1
        #vorglaetten
        varray[anzahl+1] = jacobi(varray[anzahl+1],vec(A[anzahl+1:anzahl+1, 1:end]),farray[anzahl+1],iterationen)
        #residuum
        r = farray[anzahl+1]-matrixmultiplikation(A[anzahl+1:anzahl+1, 1:end],varray[anzahl+1])
        #restriktion
        farray[anzahl] = restriktion(r)
    end

    v=zeros(1,1)
    v[1,1] = farray[anzahl][1,1]/A[anzahl,3]
    varray[anzahl] = v
    
    while(anzahl<anz)
        anzahl=anzahl+1
        #prolongation
        varray[anzahl]+= prolongation(varray[anzahl-1])
        # nachglaetten
        varray[anzahl]=jacobi(varray[anzahl],vec(A[anzahl:anzahl, 1:end]),farray[anzahl],iterationen)
    end
    return varray[anzahl]
end

function mehrgitterrbgs(A::Array{Int64,2},varray::Array{Array{Float64,2},1}, farray::Array{Array{Float64,2},1},iterationen::Int64, anzahl::Int64)
    #Rekursiver Anker
    if(anzahl==1)
        varray[anzahl][1,1]=farray[anzahl][1,1]/A[anzahl,3]
        return 0
    end
    anzahl=anzahl-1
    
    #Vorglaetten
    v=redblackgs(varray[anzahl+1],vec(A[anzahl+1:anzahl+1, 1:end]),farray[anzahl+1],iterationen)
    varray[anzahl+1]=v
    
    #residuum
    r=farray[anzahl+1]-matrixmultiplikation(A[anzahl+1:anzahl+1, 1:end],varray[anzahl+1])
    
    #restriktion
    f= restriktion(r)
    farray[anzahl]=f
    
    #rekursion
    mehrgitter(A,varray, farray,iterationen,anzahl)
    
    #prolongation
    varray[anzahl]=prolongation(varray[anzahl])
    varray[anzahl+1]+=varray[anzahl]
    #Nachglaetten
    varray[anzahl+1]=redblackgs(varray[anzahl+1],vec(A[anzahl+1:anzahl+1, 1:end]),farray[anzahl+1],iterationen)
    
end

function mehrgitter_iterativrbgs(A::Array{Int64,2},varray::Array{Array{Float64,2},1}, farray::Array{Array{Float64,2},1},iterationen::Int64, anz::Int64)
    anzahl=anz
    
    while(anzahl>1)
        anzahl=anzahl-1
        #vorglaetten
        varray[anzahl+1] = redblackgs(varray[anzahl+1],vec(A[anzahl+1:anzahl+1, 1:end]),farray[anzahl+1],iterationen)
        #residuum
        r = farray[anzahl+1]-matrixmultiplikation(A[anzahl+1:anzahl+1, 1:end],varray[anzahl+1])
        #restriktion
        farray[anzahl] = restriktion(r)
    end

    v=zeros(1,1)
    v[1,1] = farray[anzahl][1,1]/A[anzahl,3]
    varray[anzahl] = v
    
    while(anzahl<anz)
        anzahl=anzahl+1
        #prolongation
        varray[anzahl]+= prolongation(varray[anzahl-1])
        # nachglaetten
        varray[anzahl]=redblackgs(varray[anzahl],vec(A[anzahl:anzahl, 1:end]),farray[anzahl],iterationen)
    end
    return varray[anzahl]
end

#Baut auf der rekursiven Version auf mit Jacobi Glaetter
function mehrgittermethode_rekursiv(A::Array{Int64,2},varray::Array{Array{Float64,2},1}, farray::Array{Array{Float64,2},1},iterationen::Int64, anzahl::Int64,toleranz::Float64)
    i=0
    #Residuum bestimmen --> Ist Abbruchbedingung
    while(toleranz<=maximum(abs(farray[anzahl]-matrixmultiplikation(A[anzahl:anzahl, 1:end],varray[anzahl]))))
        #Alle eintraege ausser die jeweils letzten muessen genullt werden
        for j in 1:anzahl-1
            v=zeros((2^j)-1,(2^j)-1)
            varray[j]=v
            farray[j]=v
        end
        
        #Iterationsaufruf
        mehrgitter(A,varray, farray,iterationen, anzahl)
        i=i+1
    end
    return i
end

#Baut auf der iterativen Version auf mit Jacobi Glaetter
function mehrgittermethode_iterativ(A::Array{Int64,2},varray::Array{Array{Float64,2},1}, farray::Array{Array{Float64,2},1},iterationen::Int64, anzahl::Int64,toleranz::Float64)
    i=0
    #Residuum bestimmen --> Ist Abbruchbedingung
    while(toleranz<=maximum(abs(farray[anzahl]-matrixmultiplikation(A[anzahl:anzahl, 1:end],varray[anzahl]))))
        for j in 1:anzahl-1
            v=zeros((2^j)-1,(2^j)-1)
            varray[j]=v
            farray[j]=v
        end
    
        #Iterationsaufruf
        varray[anzahl]=mehrgitter_iterativ(A,varray, farray,iterationen, anzahl)
        i=i+1
    end
    return i
end

#Baut auf der rekursiven Version auf mit Jacobi Glaetter
function mehrgittermethode_rekursivrbgs(A::Array{Int64,2},varray::Array{Array{Float64,2},1}, farray::Array{Array{Float64,2},1},iterationen::Int64, anzahl::Int64,toleranz::Float64)
    i=0
    #Residuum bestimmen --> Ist Abbruchbedingung
    while(toleranz<=maximum(abs(farray[anzahl]-matrixmultiplikation(A[anzahl:anzahl, 1:end],varray[anzahl]))))
        #Alle eintraege ausser die jeweils letzten muessen genullt werden
        for j in 1:anzahl-1
            v=zeros((2^j)-1,(2^j)-1)
            varray[j]=v
            farray[j]=v
        end
        
        #Iterationsaufruf
        mehrgitterrbgs(A,varray, farray,iterationen, anzahl)
        i=i+1
    end
    return i
end

#Baut auf der iterativen Version auf mit Jacobi Glaetter
function mehrgittermethode_iterativrbgs(A::Array{Int64,2},varray::Array{Array{Float64,2},1}, farray::Array{Array{Float64,2},1},iterationen::Int64, anzahl::Int64,toleranz::Float64)
    i=0
    #Residuum bestimmen --> Ist Abbruchbedingung
    while(toleranz<=maximum(abs(farray[anzahl]-matrixmultiplikation(A[anzahl:anzahl, 1:end],varray[anzahl]))))
        for j in 1:anzahl-1
            v=zeros((2^j)-1,(2^j)-1)
            varray[j]=v
            farray[j]=v
        end
    
        #Iterationsaufruf
        varray[anzahl]=mehrgitter_iterativrbgs(A,varray, farray,iterationen, anzahl)
        i=i+1
    end
    return i
end

function loeseoo(farray::Array{Hmatrix,1}, A::Array{Stencil,1},level::Int64)
    v=Hmatrix(zeros(1,1))
    v.element[1,1] = farray[level].element[1][1]/A[level].element[3]
    return v
end
#erzeugt einen 2Dimensionalen Array --> Fouriermode
function fourier(N::Int64,k::Int64)
    f(x,y)=sin(k*pi*x)*sin(k*pi*y)
    matrix=zeros(N,N)
    for i in 1:N
        for j in 1:N
            matrix[i,j]=f(i/(N+1),j/(N+1)) 
        end
    end
    return matrix
end

# Matrixmultiplikation mit einem Stencil im zweidimensionalen Bereich
function matrixmultiplikationoo(A::Stencil,un::Hmatrix)
    u=fuege_Rand_hinzuoo(un)
    gr=size(u.element)[1]
    matrix=Hmatrix(zeros(gr-2,gr-2))
    for i in 2:gr-1
        for j in 2:gr-1
            matrix.element[i-1,j-1]=A.element[1]*u.element[i-1,j]+A.element[2]*u.element[i,j-1]+A.element[3]*u.element[i,j]+A.element[4]*u.element[i,j+1]+A.element[5]*u.element[i+1,j]
        end
    end
    return matrix
end

#Ich arbeite mit einer geraenderten Matrix bei diesem Verfahren,
#dies macht die Implementierung leichter.
function fuege_Rand_hinzuoo(u::Hmatrix)
    N=size(u.element)[1]
    matrix=Hmatrix(zeros(N+2,N+2))
    matrix.element[2:N+1,2:N+1] = u.element
    return matrix
end

#2-Dimensionaler Jacobi Algorithmus
#Schau in die Doku wenn du nicht weisst wies funktioniert!
function jacobioo(x::Hmatrix,a::Stencil,b::Hmatrix,anzahl::Int64)
    w=2/3
    x=fuege_Rand_hinzuoo(x)
    b=fuege_Rand_hinzuoo(b)
    n=size(x.element)[1]
    z=0
    while(z<anzahl)
        xalt=Hmatrix(copy(x.element))
        for i in 2:n-1
            for j in 2:n-1
                x.element[i,j]= (1-w)*xalt.element[i,j]+((w/a.element[3])*(b.element[i,j]-(a.element[1]*xalt.element[i-1,j]+a.element[2]*xalt.element[i,j-1]+a.element[4]*xalt.element[i,j+1]+a.element[5]*xalt.element[i+1,j])))
            end
        end
        z=z+1
    end
    x.element=x.element[2:end-1, 2:end-1]
    return x
end

function redblackgsoo(x::Hmatrix,a::Stencil,b::Hmatrix,anzahl::Int64)
    x=fuege_Rand_hinzuoo(x)
    b=fuege_Rand_hinzuoo(b)
    n=size(x.element)[1]
    z=0
    while(z<anzahl)
        #Erste Schleife
        for i in 2:(n-1)
            for j in collect(2+(i%2):2:n-1)
                x.element[i,j]= ((1/a.element[3])*(b.element[i,j]-(a.element[1]*x.element[i-1,j]+a.element[2]*x.element[i,j-1]+a.element[4]*x.element[i,j+1]+a.element[5]*x.element[i+1,j])))
            end
        end
        #Zweite Schleife
        for i in 2:(n-1)
            for j in collect(3-(i%2):2:n-1)
                x.element[i,j]= ((1/a.element[3])*(b.element[i,j]-(a.element[1]*x.element[i-1,j]+a.element[2]*x.element[i,j-1]+a.element[4]*x.element[i,j+1]+a.element[5]*x.element[i+1,j])))
            end
        end
        z=z+1
    end
    x.element=x.element[2:end-1, 2:end-1]
    return x
end

#Auch hier arbeite ich mit einer geraenderten Matrix
#Die Matrix wird normal restringiert
function restriktionoo(x::Hmatrix)
    groesse = size(x.element)[1]
    matrix=Hmatrix(zeros(convert(Int,(groesse-1)/2),convert(Int,(groesse-1)/2)))
    x=fuege_Rand_hinzuoo(x)
    k=3
    for i in 1:convert(Int,(groesse-1)/2)
        l=3
        for j in 1:convert(Int,(groesse-1)/2)
            matrix.element[i,j]=0.5*x.element[k,l]+0.125*x.element[k-1,l]+0.125*x.element[k+1,l]+0.125*x.element[k,l-1]+0.125*x.element[k,l+1]
            l=l+2
        end
        k=k+2
    end
    return matrix 
end

#Die feinere Matrix wird in die groebere Matrix gesetzt
function einsetzenoo(matrix::Hmatrix)
    groesse= size(matrix.element)[1]
    neue=Hmatrix(zeros(groesse*2+3, groesse*2+3))
    k=1
    for i in range(3,2,groesse)
        l=1
        for j in range(3,2,groesse)
            neue.element[i,j]=matrix.element[k,l]
            l=l+1
        end
        k=k+1
     end
    return neue
end

#Die Matrix wird beim einsetzen geraendert und es erfolgt eine Prolongation
function prolongationoo(xi::Hmatrix)
    x=einsetzenoo(xi)
    groesse= size(x.element)[1]
    #1. Zeilen in denen schon Werte stehen werden bearbeitet
    for i in range(3,2,convert(Int,(groesse-3)/2))
        for j in range(2,2,convert(Int,(groesse-3)/2+1))
            x.element[i,j]=0.5*x.element[i, j-1]+0.5*x.element[i, j+1]
        end
    end
    #2. Zeilen in denen noch keine Werte stehen werden bearbeitet
    for i in range(2,2,convert(Int,(groesse-3)/2+1))
        for j in range(2,1,groesse-2)
            x.element[i,j]=0.5*x.element[i-1, j]+0.5*x.element[i+1, j]
        end
    end
    xi.element=x.element[2:end-1, 2:end-1]
    return xi
end

function mehrgitteroo(A::Array{Stencil},varray::Array{Hmatrix}, farray::Array{Hmatrix},iterationen::Int64, level::Int64)
    #Rekursiver Anker
    if(level==1)
        varray[level]=loeseoo(farray,A,level)
        return 0
    end
    level=level-1
    #Vorglaetten
    varray[level+1]=jacobioo(varray[level+1],A[level+1],farray[level+1],iterationen)
    
    #residuum
    r=farray[level+1]-matrixmultiplikationoo(A[level+1],varray[level+1])
    
    #restriktion
    f= restriktionoo(r)
    farray[level]=f
    
    #rekursion
    mehrgitteroo(A,varray, farray,iterationen,level)
    
    #prolongation
    varray[level]=prolongationoo(varray[level])
    varray[level+1]=varray[level+1]+varray[level]
    #Nachglaetten
    varray[level+1]=jacobioo(varray[level+1],A[level+1],farray[level+1],iterationen)
    
end

function mehrgitter_iterativoo(A::Array{Stencil},varray::Array{Hmatrix}, farray::Array{Hmatrix},iterationen::Int64, anz::Int64)
    level=anz
    
    while(level>1)
        level=level-1
        #vorglaetten
        varray[level+1]= jacobioo(varray[level+1],A[level+1],farray[level+1],iterationen)
        #residuum
        r = farray[level+1]-matrixmultiplikationoo(A[level+1],varray[level+1])
        #restriktion
        farray[level] = restriktionoo(r)
    end
    
    varray[level]= loeseoo(farray, A,level)
    
    while(level<anz)
        level=level+1
        #prolongation
        varray[level]= varray[level]+prolongationoo(varray[level-1])
        # nachglaetten
        varray[level]=jacobioo(varray[level],A[level],farray[level],iterationen)
    end
    return varray[level]
end

function mehrgitteroorbgs(A::Array{Stencil},varray::Array{Hmatrix}, farray::Array{Hmatrix},iterationen::Int64, level::Int64)
    #Rekursiver Anker
    if(level==1)
        varray[level]=loeseoo(farray,A,level)
        return 0
    end
    level=level-1
    #Vorglaetten
    varray[level+1]=redblackgsoo(varray[level+1],A[level+1],farray[level+1],iterationen)
    
    #residuum
    r=farray[level+1]-matrixmultiplikationoo(A[level+1],varray[level+1])
    
    #restriktion
    f= restriktionoo(r)
    farray[level]=f
    
    #rekursion
    mehrgitteroo(A,varray, farray,iterationen,level)
    
    #prolongation
    varray[level]=prolongationoo(varray[level])
    varray[level+1]=varray[level+1]+varray[level]
    #Nachglaetten
    varray[level+1]=redblackgsoo(varray[level+1],A[level+1],farray[level+1],iterationen)
    
end

function mehrgitter_iterativoorbgs(A::Array{Stencil},varray::Array{Hmatrix}, farray::Array{Hmatrix},iterationen::Int64, anz::Int64)
    level=anz
    
    while(level>1)
        level=level-1
        #vorglaetten
        varray[level+1]= redblackgsoo(varray[level+1],A[level+1],farray[level+1],iterationen)
        #residuum
        r = farray[level+1]-matrixmultiplikationoo(A[level+1],varray[level+1])
        #restriktion
        farray[level] = restriktionoo(r)
    end
    
    varray[level]= loeseoo(farray, A,level)
    
    while(level<anz)
        level=level+1
        #prolongation
        varray[level]= varray[level]+prolongationoo(varray[level-1])
        # nachglaetten
        varray[level]=redblackgsoo(varray[level],A[level],farray[level],iterationen)
    end
    return varray[level]
end

#Baut auf der rekursiven Version auf mit Jacobi Glaetter
function mehrgittermethode_rekursivoo(A::Array{Stencil},varray::Array{Hmatrix}, farray::Array{Hmatrix},iterationen::Int64, level::Int64,toleranz::Float64)
    i=0
    #Residuum bestimmen --> Ist Abbruchbedingung
    while(toleranz<=maximum(abs(farray[level]-matrixmultiplikationoo(A[level],varray[level]))))
        #Alle eintraege ausser die jeweils letzten muessen genullt werden
        for j in 1:level-1
            v=Hmatrix(zeros((2^j)-1,(2^j)-1))
            varray[j]=v
            farray[j]=v
        end
        
        #Iterationsaufruf
        mehrgitteroo(A,varray, farray,iterationen, level)
        i=i+1
    end
    return i
end

#Baut auf der iterativen Version auf mit Jacobi Glaetter
function mehrgittermethode_iterativoorbgs(A::Array{Stencil},varray::Array{Hmatrix}, farray::Array{Hmatrix},iterationen::Int64, level::Int64,toleranz::Float64)
    i=0
    #Residuum bestimmen --> Ist Abbruchbedingung
    while(toleranz<=maximum(abs(farray[level]-matrixmultiplikationoo(A[level],varray[level]))))
        for j in 1:level-1
            v=Hmatrix(zeros((2^j)-1,(2^j)-1))
            varray[j]=v
            farray[j]=v
        end
    
        #Iterationsaufruf
        varray[level]=mehrgitter_iterativoorbgs(A,varray, farray,iterationen, level)
        i=i+1
    end
    return i
end

#Baut auf der rekursiven Version auf mit Jacobi Glaetter
function mehrgittermethode_rekursivoorbgs(A::Array{Stencil},varray::Array{Hmatrix}, farray::Array{Hmatrix},iterationen::Int64, level::Int64,toleranz::Float64)
    i=0
    #Residuum bestimmen --> Ist Abbruchbedingung
    while(toleranz<=maximum(abs(farray[level]-matrixmultiplikationoo(A[level],varray[level]))))
        #Alle eintraege ausser die jeweils letzten muessen genullt werden
        for j in 1:level-1
            v=Hmatrix(zeros((2^j)-1,(2^j)-1))
            varray[j]=v
            farray[j]=v
        end
        
        #Iterationsaufruf
        mehrgitteroorbgs(A,varray, farray,iterationen, level)
        i=i+1
    end
    return i
end

#Baut auf der iterativen Version auf mit Jacobi Glaetter
function mehrgittermethode_iterativoo(A::Array{Stencil},varray::Array{Hmatrix}, farray::Array{Hmatrix},iterationen::Int64, level::Int64,toleranz::Float64)
    i=0
    #Residuum bestimmen --> Ist Abbruchbedingung
    while(toleranz<=maximum(abs(farray[level]-matrixmultiplikationoo(A[level],varray[level]))))
        for j in 1:level-1
            v=Hmatrix(zeros((2^j)-1,(2^j)-1))
            varray[j]=v
            farray[j]=v
        end
    
        #Iterationsaufruf
        varray[level]=mehrgitter_iterativoo(A,varray, farray,iterationen, level)
        i=i+1
    end
    return i
end

n=31
toleranz=0.1

ausgabe = open("/home/zam/heisters/projects/iterative_loeser/Visualisierung/Visualisierung_Daten/ergebnis_mehrgitter_zeit.txt", "w")

while(n<10000)
    println("Anzahl:",n)
    #mache AArray
    anzahl=convert(Int64,log2(n+1))
    aarray =Array{Int64,2}(anzahl,5)
    array = [1,1,-4,1,1]
    for i in 1:5
        for j in 1:anzahl
            aarray[j,i]=-(2^j)^2*array[i]
        end
    end
    schritte = 0
    zeit = zeros(25)
    for w in 1:25
        #Mache f und v
        srand(12345)
        varray=Array{Array{Float64,2}}(anzahl)
        varray[anzahl]=abs(sin(rand(n,n)))
        farray=Array{Array{Float64,2}}(anzahl)
        farray[anzahl]=zeros(n,n)
        t=time()
        schritte=mehrgittermethode_rekursiv(aarray,varray, farray,2, anzahl, toleranz)
        zeit[w]=time()-t
    end
    println(sum(zeit)/size(zeit)[1])
    write(ausgabe, string(n),";",string(sum(zeit)/size(zeit)[1]),";",string(schritte),";rek;\n")
    
    for w in 1:25
        srand(12345)
        varray=Array{Array{Float64,2}}(anzahl)
        varray[anzahl]=abs(sin(rand(n,n)))
        farray=Array{Array{Float64,2}}(anzahl)
        farray[anzahl]=zeros(n,n)
        t=time()
        schritte= mehrgittermethode_iterativ(aarray,varray, farray,2, anzahl, toleranz)
        zeit[w]=time()-t
    end
    println(sum(zeit)/size(zeit)[1])
    write(ausgabe, string(n),";",string(sum(zeit)/size(zeit)[1]),";",string(schritte),";itv;\n")
    n=2*(n+1)-1
end
close(ausgabe)

ausgabe = open("/home/zam/heisters/projects/iterative_loeser/Visualisierung/Visualisierung_Daten/ergebnis_mehrgitter_objektorientiert_zeit.txt", "w")
n=31
while(n<10000)
    println(n)
    
    #Mache A Array
    level=convert(Int64,log2(n+1))
    aarray =Array{Stencil}(level)
    for i in 1:level
        aarray[i]=Stencil(-(2^i)^2*[1,1,-4,1,1])
    end
   
    schritte = 0
    
    zeit = zeros(25)
    for w in 1:25
        srand(12345)
        f=Hmatrix(zeros(n,n))
        v=Hmatrix(abs(sin(rand(n,n))))
        varray=Array{Hmatrix}(level)
        varray[level]=v
        farray=Array{Hmatrix}(level)
        farray[level]=f
        t=time()
        #Mache f und v
        schritte=mehrgittermethode_rekursivoo(aarray,varray, farray,2, level, toleranz)
        zeit[w]=time()-t
    end
    println("rekursiv:", sum(zeit)/size(zeit)[1])
    write(ausgabe, string(n),";",string(sum(zeit)/size(zeit)[1]),";",string(schritte),";rek;\n")
    
    for w in 1:25
        srand(12345)
        f=Hmatrix(zeros(n,n))
        v=Hmatrix(abs(sin(rand(n,n))))
        varray=Array{Hmatrix}(level)
        varray[level]=v
        farray=Array{Hmatrix}(level)
        farray[level]=f
        t=time()
        schritte= mehrgittermethode_iterativoo(aarray,varray, farray,2, level, toleranz)
        zeit[w]=time()-t
    end
    println("iterativ:",sum(zeit)/size(zeit)[1])
    write(ausgabe, string(n),";",string(sum(zeit)/size(zeit)[1]),";",string(schritte),";itv;\n")
    n=2*(n+1)-1
end
close(ausgabe)

ausgabe = open("/home/zam/heisters/projects/iterative_loeser/Visualisierung/Visualisierung_Daten/ergebnis_mehrgitter_zeit_rbgs.txt", "w")
n=31
while(n<10000)
    println("Anzahl:",n)
    #mache AArray
    anzahl=convert(Int64,log2(n+1))
    aarray =Array{Int64,2}(anzahl,5)
    array = [1,1,-4,1,1]
    for i in 1:5
        for j in 1:anzahl
            aarray[j,i]=-(2^j)^2*array[i]
        end
    end
    schritte = 0
    zeit = zeros(25)
    for w in 1:25
        #Mache f und v
        srand(12345)
        varray=Array{Array{Float64,2}}(anzahl)
        varray[anzahl]=abs(sin(rand(n,n)))
        farray=Array{Array{Float64,2}}(anzahl)
        farray[anzahl]=zeros(n,n)
        t=time()
        schritte=mehrgittermethode_rekursiv(aarray,varray, farray,2, anzahl, toleranz)
        zeit[w]=time()-t
    end
    println(sum(zeit)/size(zeit)[1])
    write(ausgabe, string(n),";",string(sum(zeit)/size(zeit)[1]),";",string(schritte),";rek;\n")
    
    for w in 1:25
        srand(12345)
        varray=Array{Array{Float64,2}}(anzahl)
        varray[anzahl]=abs(sin(rand(n,n)))
        farray=Array{Array{Float64,2}}(anzahl)
        farray[anzahl]=zeros(n,n)
        t=time()
        schritte= mehrgittermethode_iterativ(aarray,varray, farray,2, anzahl, toleranz)
        zeit[w]=time()-t
    end
    println(sum(zeit)/size(zeit)[1])
    write(ausgabe, string(n),";",string(sum(zeit)/size(zeit)[1]),";",string(schritte),";itv;\n")
    n=2*(n+1)-1
end
close(ausgabe)

ausgabe = open("/home/zam/heisters/projects/iterative_loeser/Visualisierung/Visualisierung_Daten/ergebnis_mehrgitter_objektorientiert_zeit_rbgs.txt", "w")
n=31
while(n<10000)
    println(n)
    
    #Mache A Array
    level=convert(Int64,log2(n+1))
    aarray =Array{Stencil}(level)
    for i in 1:level
        aarray[i]=Stencil(-(2^i)^2*[1,1,-4,1,1])
    end
   
    schritte = 0
    
    zeit = zeros(25)
    for w in 1:25
        srand(12345)
        f=Hmatrix(zeros(n,n))
        v=Hmatrix(abs(sin(rand(n,n))))
        varray=Array{Hmatrix}(level)
        varray[level]=v
        farray=Array{Hmatrix}(level)
        farray[level]=f
        t=time()
        #Mache f und v
        schritte=mehrgittermethode_rekursivoo(aarray,varray, farray,2, level, toleranz)
        zeit[w]=time()-t
    end
    println("rekursiv:", sum(zeit)/size(zeit)[1])
    write(ausgabe, string(n),";",string(sum(zeit)/size(zeit)[1]),";",string(schritte),";rek;\n")
    
    for w in 1:25
        srand(12345)
        f=Hmatrix(zeros(n,n))
        v=Hmatrix(abs(sin(rand(n,n))))
        varray=Array{Hmatrix}(level)
        varray[level]=v
        farray=Array{Hmatrix}(level)
        farray[level]=f
        t=time()
        schritte= mehrgittermethode_iterativoo(aarray,varray, farray,2, level, toleranz)
        zeit[w]=time()-t
    end
    println("iterativ:",sum(zeit)/size(zeit)[1])
    write(ausgabe, string(n),";",string(sum(zeit)/size(zeit)[1]),";",string(schritte),";itv;\n")
    n=2*(n+1)-1
end
close(ausgabe)

