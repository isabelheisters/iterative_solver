using MPI

#erzeugt einen 2Dimensionalen Array --> Fouriermode
function fourier(f::Function, N::Int64,k::Int64, globx::Int64, globy::Int64, anzx::Int64, anzy::Int64)
    matrix=zeros(anzx,anzy)
    l=1
    for i in globx:globx+anzx-1
        d=1
        for j in globy:globy+anzy-1
            matrix[l,d]=f(i/(N+1),j/(N+1)) 
            d=d+1
        end
        l=l+1
    end
    return matrix
end

# Matrixmultiplikation mit einem Stencil im zweidimensionalen Bereich
function matrixmultiplikation(A::Array{Int64,1},u::Array{Float64,2})
    gr1=size(u)[1]
    gr2=size(u)[2]
    matrix=zeros(gr1-2,gr2-2)
    for i in 2:gr1-1
        for j in 2:gr2-1
            matrix[i-1,j-1]=A[1]*u[i-1,j]+A[2]*u[i,j-1]+A[3]*u[i,j]+A[4]*u[i,j+1]+A[5]*u[i+1,j]
        end
    end
    return matrix
end

#Eine geraenderte Matrix wird uebergeben. Von den nebenliegenden Prozessoren werden die raender gefuellt
function halo_austausch!(id::Int64, size::Int64, x::Array{Float64,2})
    #bei einem prozessor garnichts machen
    if size!=1
        
        #alle sachen die ausgetauscht werden muessen
        oben=x[2:2,2:anzy+1]
        unten=x[anzx+1:anzx+1, 2:anzy+1]
        links=x[2:anzx+1]
        rechts=x[2:anzx+1, anzy+1:anzy+1]
        obenbek=zeros(oben)
        untenbek=zeros(unten)
        linksbek=zeros(links)
        rechtsbek=zeros(rechts)
        
        #Als erstes wird spaltenweise getauscht nach rechts und nach links, zuerst senden alle geraden spalten, dann #die ungeraden
        if id%sqrt(size)%2==0
            #in der ersten spalte wird nur nach rechts getauscht
            if id%sqrt(size)==0 
                MPI.Send(rechts,id+1,0,MPI.COMM_WORLD)
                MPI.Recv!(rechtsbek,id+1,0,MPI.COMM_WORLD)
                x[2:anzx+1,anzy+2:anzy+2]=rechtsbek
            #in der letzten spalte wird nur nach links getauscht
            elseif id%sqrt(size)==sqrt(size)-1
                MPI.Send(links,id-1,0,MPI.COMM_WORLD)
                MPI.Recv!(linksbek,id-1,0,MPI.COMM_WORLD)
                x[2:anzx+1,1:1]=linksbek
            else
                MPI.Send(rechts,id+1,0,MPI.COMM_WORLD)
                MPI.Send(links,id-1,0,MPI.COMM_WORLD)
                MPI.Recv!(rechtsbek,id+1,0,MPI.COMM_WORLD)
                x[2:anzx+1,anzy+2:anzy+2]=rechtsbek
                MPI.Recv!(linksbek,id-1,0,MPI.COMM_WORLD)
                x[2:anzx+1,1:1]=linksbek
            end
        else
            #in der letzten spalte wird nur nach links getauscht
            if id%sqrt(size)==sqrt(size)-1
                MPI.Recv!(linksbek,id-1,0,MPI.COMM_WORLD)
                x[2:anzx+1,1:1]=linksbek
                MPI.Send(links,id-1,0,MPI.COMM_WORLD)
            else
                MPI.Recv!(rechtsbek,id+1,0,MPI.COMM_WORLD)
                x[2:anzx+1,anzy+2:anzy+2]=rechtsbek
                MPI.Recv!(linksbek,id-1,0,MPI.COMM_WORLD)
                x[2:anzx+1,1:1]=linksbek
                MPI.Send(rechts,id+1,0,MPI.COMM_WORLD)
                MPI.Send(links,id-1,0,MPI.COMM_WORLD)
            end
        end
        
        #jetzt wird spaltenweise getauscht, erst die geraden spalten senden dann die ungeraden
        if floor(id/sqrt(size))%2==0
            #in der ersten zeile wird nur nach unten getauscht
            if floor(id/sqrt(size))==0 
                MPI.Send(unten,convert(Int32,id+sqrt(size)),0,MPI.COMM_WORLD)
                MPI.Recv!(untenbek,convert(Int32,id+sqrt(size)),0,MPI.COMM_WORLD)
                x[anzx+2:anzx+2,2:anzy+1]=untenbek
            #in der letzten zeile wird nur nach oben getauscht
            elseif floor(id/sqrt(size))==sqrt(size)-1
                MPI.Send(oben,convert(Int32,id-sqrt(size)),0,MPI.COMM_WORLD)
                MPI.Recv!(obenbek,convert(Int32,id-sqrt(size)),0,MPI.COMM_WORLD)
                x[1:1,2:anzy+1]=obenbek
            else
                MPI.Send(unten,convert(Int32,id+sqrt(size)),0,MPI.COMM_WORLD)
                MPI.Send(oben,convert(Int32,id-sqrt(size)),0,MPI.COMM_WORLD)
                MPI.Recv!(untenbek,convert(Int32,id+sqrt(size)),0,MPI.COMM_WORLD)
                x[anzx+2:anzx+2,2:anzy+1]=untenbek
                MPI.Recv!(obenbek,convert(Int32,id-sqrt(size)),0,MPI.COMM_WORLD)
                x[1:1,2:anzy+1]=obenbek
            end
        else
            #in der letzten zeile wird nur nach oben getauscht
            if floor(id/sqrt(size))==sqrt(size)-1
                MPI.Recv!(obenbek,convert(Int32,id-sqrt(size)),0,MPI.COMM_WORLD)
                x[1:1,2:anzy+1]=obenbek
                MPI.Send(oben,convert(Int32,id-sqrt(size)),0,MPI.COMM_WORLD)
            else
                MPI.Recv!(untenbek,convert(Int32,id+sqrt(size)),0,MPI.COMM_WORLD)
                x[anzx+2:anzx+2,2:anzy+1]=untenbek
                MPI.Recv!(obenbek,convert(Int32,id-sqrt(size)),0,MPI.COMM_WORLD)
                x[1:1,2:anzy+1]=obenbek
                MPI.Send(unten,convert(Int32,id+sqrt(size)),0,MPI.COMM_WORLD)
                MPI.Send(oben,convert(Int32,id-sqrt(size)),0,MPI.COMM_WORLD)
            end
        end
    end 
end

#Ich arbeite mit einer geraenderten Matrix bei diesem Verfahren,
#dies macht die Implementierung leichter.
function fuege_Rand_hinzu(u::Array{Float64,2})
    N1=size(u)[1]
    N2=size(u)[2]
    matrix=zeros(N1+2,N2+2)
    matrix[2:N1+1,2:N2+1] = u
    return matrix
end

function main(n::Int64, k::Int64)
   
    #es wird herausgefynden wie viele prozessore und welche id diese haben
    id = MPI.Comm_rank(MPI.COMM_WORLD)
    size = MPI.Comm_size(MPI.COMM_WORLD)
    
    globx=1
    globy=1
    anzx=0
    anzy=0
    
    #die prozessorenanzahl darf nur eine quadratzahl sein
    if size==0 && convert(Int32, sqrt(size))^2!=size 
        if size!=1
            println(size!=0 && convert(Int32, sqrt(size))^2==size)
            return 0
        end
    end
    
    #Die globalen x werte werden bestimmt und die matrixgroesse und das bedingt pro prozessor
    globy=convert(Int64,ceil((n+1)/sqrt(size))*(id%sqrt(size))+1)
    if id%sqrt(size)==sqrt(size)-1
        anzy=convert(Int64,n-globy+1)
    else
        anzy=convert(Int64,ceil((n+1)/sqrt(size)))
    end
    
    globx=convert(Int64,floor(id/sqrt(size))*ceil((n+1)/sqrt(size))+1)
    
    if id>=size-sqrt(size)
        anzx=convert(Int64,n-globx+1)
    else
        anzx=convert(Int64, ceil((n+1)/sqrt(size)))
    end
    
    #Funktionen und A und x werden erzeugt
    f(x,y)=sin(k*pi*x)*sin(k*pi*y)
    f2(x,y)=-2*(k*k*pi*pi)*sin(k*pi*x)*sin(k*pi*y)  
    A=[((n+1)*(n+1)), ((n+1)*(n+1)), -4*((n+1)*(n+1)), ((n+1)*(n+1)), ((n+1)*(n+1))]
    x=fourier(f,n,k,globx,globy,anzx,anzy)

    #Rand hinzufuegen
    x=fuege_Rand_hinzu(x)
    
    #Haloaustausch
    halo_austausch!(id,size,x)
    
    #Matrixmultiplikation
    x=matrixmultiplikation(A,x)
    
    #Maximum bestimmen
    max=maximum(abs(fourier(f2,n,k,globx,globy,anzx,anzy)-x))
    y=MPI.Reduce(max, MPI.MAX, 0, MPI.COMM_WORLD)
    if(id==0)
        println("Error: ",y)
    end
    return y
end

MPI.Init()
id = MPI.Comm_rank(MPI.COMM_WORLD)
n=7
k=1
if(id==0)
   datei=open("Daten/matrixmultiplikation_zweidimensional_4_threads.txt","w")
end
while n<10000
    error=main(n,k)
    if(id==0)
        write(datei,string(n),";",string(k),";"string(error), ";\n")
    end
    n=2*(n+1)-1
end
if(id==0)
    close(datei)
end
MPI.Finalize()
    
