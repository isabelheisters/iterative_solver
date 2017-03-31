using MPI

function redblackgs(id::Int64, num_procs::Int64, x::Array{Float64,2},a::Array{Int64,1},b::Array{Float64,2},anzahl::Int64)
    x=fuege_Rand_hinzu(x)
    b=fuege_Rand_hinzu(b)
    gr1=size(x)[1]
    gr2=size(x)[2]
    z=0
    while(z<anzahl)
        #Erste Schleife
        halo_austausch!(id, num_procs,x,gr1-2, gr2-2)
        if num_procs%2==0
            for i in 2:(gr1-1)
                for j in collect(2+(i%2):2:gr2-1)
                    x[i,j]= ((1/a[3])*(b[i,j]-(a[1]*x[i-1,j]+a[2]*x[i,j-1]+a[4]*x[i,j+1]+a[5]*x[i+1,j])))
                end
            end
            #Zweite Schleife
            halo_austausch!(id, num_procs,x,gr1-2, gr2-2)
            for i in 2:(gr1-1)
                for j in collect(3-(i%2):2:gr2-1)
                    x[i,j]= ((1/a[3])*(b[i,j]-(a[1]*x[i-1,j]+a[2]*x[i,j-1]+a[4]*x[i,j+1]+a[5]*x[i+1,j])))
                end
            end
        else
            for i in 2:(gr1-1)
                for j in collect(2+(id%2):2:gr2-1)
                    x[i,j]= ((1/a[3])*(b[i,j]-(a[1]*x[i-1,j]+a[2]*x[i,j-1]+a[4]*x[i,j+1]+a[5]*x[i+1,j])))
                end
            end
            #Zweite Schleife
            halo_austausch!(id, num_procs,x,gr1-2, gr2-2)
            for i in 2:(gr1-1)
                for j in collect(3-(id%2):2:gr2-1)
                    x[i,j]= ((1/a[3])*(b[i,j]-(a[1]*x[i-1,j]+a[2]*x[i,j-1]+a[4]*x[i,j+1]+a[5]*x[i+1,j])))
                end
            end
        end
        z=z+1
    end
    x=x[2:end-1, 2:end-1]
    return x
end

function residuumsbedingung(A::Array{Int64,1}, x::Array{Float64,2}, b::Array{Float64,2})
    y=[0.0]
    max=maximum(abs(b-matrixmultiplikation(A,x)))
    MPI.Allreduce!([max],y,MPI.MAX, MPI.COMM_WORLD)
    return y[1]
end

function redblackgs_mit_residuum(id::Int64, num_procs::Int64, x::Array{Float64,2},a::Array{Int64,1},b::Array{Float64,2},residuum::Float64)
    x=fuege_Rand_hinzu(x)
    b=fuege_Rand_hinzu(b)
    gr1=size(x)[1]
    gr2=size(x)[2]
    halo_austausch!(id, num_procs,x,gr1-2, gr2-2)
    while(residuum<=residuumsbedingung(a,x,b))
        #Erste Schleife
        if num_procs%2==0
            for i in 2:(gr1-1)
                for j in collect(2+(i%2):2:gr2-1)
                    x[i,j]= ((1/a[3])*(b[i,j]-(a[1]*x[i-1,j]+a[2]*x[i,j-1]+a[4]*x[i,j+1]+a[5]*x[i+1,j])))
                end
            end
            #Zweite Schleife
            halo_austausch!(id, num_procs,x,gr1-2, gr2-2)
            for i in 2:(gr1-1)
                for j in collect(3-(i%2):2:gr2-1)
                    x[i,j]= ((1/a[3])*(b[i,j]-(a[1]*x[i-1,j]+a[2]*x[i,j-1]+a[4]*x[i,j+1]+a[5]*x[i+1,j])))
                end
            end
            halo_austausch!(id, num_procs,x,gr1-2, gr2-2)
        else
            for i in 2:(gr1-1)
                for j in collect(2+(id%2):2:gr2-1)
                    x[i,j]= ((1/a[3])*(b[i,j]-(a[1]*x[i-1,j]+a[2]*x[i,j-1]+a[4]*x[i,j+1]+a[5]*x[i+1,j])))
                end
            end
            #Zweite Schleife
            halo_austausch!(id, num_procs,x,gr1-2, gr2-2)
            for i in 2:(gr1-1)
                for j in collect(3-(id%2):2:gr2-1)
                    x[i,j]= ((1/a[3])*(b[i,j]-(a[1]*x[i-1,j]+a[2]*x[i,j-1]+a[4]*x[i,j+1]+a[5]*x[i+1,j])))
                end
            end
            halo_austausch!(id, num_procs,x,gr1-2, gr2-2)
        end
    end
    x=x[2:end-1, 2:end-1]
    return x
end

function matrixmultiplikation(A::Array{Int64,1},u::Array{Float64,2})
    gr1=size(u)[1]
    gr2=size(u)[2]
    matrix=zeros(gr1,gr2)
    for i in 2:gr1-1
        for j in 2:gr2-1
            matrix[i,j]=A[1]*u[i-1,j]+A[2]*u[i,j-1]+A[3]*u[i,j]+A[4]*u[i,j+1]+A[5]*u[i+1,j]
        end
    end
    return matrix
end

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

#Eine geraenderte Matrix wird uebergeben. Von den nebenliegenden Prozessoren werden die raender gefuellt
function halo_austausch!(id::Int64, num_procs::Int64, x::Array{Float64,2}, anzx::Int64, anzy::Int64)
    #bei einem prozessor garnichts machen
    if num_procs!=1
    
        if floor(sqrt(num_procs))^2==num_procs
            xprocs=convert(Int64, sqrt(num_procs))
            yprocs=xprocs
        else
            xprocs=sqrt(convert(Int64, num_procs/2))
            yprocs=xprocs*2
        end
        
        #alle sachen die ausgetauscht werden muessen
        oben=x[2:2,2:anzy+1]
        unten=x[anzx+1:anzx+1, 2:anzy+1]
        links=x[2:anzx+1, 2:2]
        rechts=x[2:anzx+1, anzy+1:anzy+1]
        obenbek=zeros(oben)
        untenbek=zeros(unten)
        linksbek=zeros(links)
        rechtsbek=zeros(rechts)
        
        #Als erstes wird spaltenweise getauscht nach rechts und nach links, zuerst senden alle geraden spalten, dann die ungeraden
        if num_procs!=2    
            if id%xprocs%2==0
                #in der ersten spalte wird nur nach rechts getauscht
                if id%xprocs==0 
                    MPI.Send(rechts,id+1,0,MPI.COMM_WORLD)
                    MPI.Recv!(rechtsbek,id+1,0,MPI.COMM_WORLD)
                    x[2:anzx+1,anzy+2:anzy+2]=rechtsbek
                #in der letzten spalte wird nur nach links getauscht
                elseif id%xprocs==xprocs-1
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
                if id%xprocs==xprocs-1
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
        end
        #jetzt wird zeilenweise getauscht, erst die geraden zeilen senden dann die ungeraden
        if floor(id/xprocs)%2==0
            #in der ersten zeile wird nur nach unten getauscht
            if floor(id/xprocs)==0 
                MPI.Send(unten,convert(Int64,id+xprocs),0,MPI.COMM_WORLD)
                MPI.Recv!(untenbek,convert(Int64,id+xprocs),0,MPI.COMM_WORLD)
                x[anzx+2:anzx+2,2:anzy+1]=untenbek
            #in der letzten zeile wird nur nach oben getauscht
            elseif floor(id/xprocs)==yprocs-1
                MPI.Send(oben,convert(Int64,id-xprocs),0,MPI.COMM_WORLD)
                MPI.Recv!(obenbek,convert(Int64,id-xprocs),0,MPI.COMM_WORLD)
                x[1:1,2:anzy+1]=obenbek
            else
                MPI.Send(unten,convert(Int64,id+xprocs),0,MPI.COMM_WORLD)
                MPI.Send(oben,convert(Int64,id-xprocs),0,MPI.COMM_WORLD)
                MPI.Recv!(untenbek,convert(Int64,id+xprocs),0,MPI.COMM_WORLD)
                x[anzx+2:anzx+2,2:anzy+1]=untenbek
                MPI.Recv!(obenbek,convert(Int64,id-xprocs),0,MPI.COMM_WORLD)
                x[1:1,2:anzy+1]=obenbek
            end
        else
            
            #in der letzten zeile wird nur nach oben getauscht
            if floor(id/xprocs)==yprocs-1
                MPI.Recv!(obenbek,convert(Int64,id-xprocs),0,MPI.COMM_WORLD)
                x[1:1,2:anzy+1]=obenbek
                MPI.Send(oben,convert(Int64,id-xprocs),0,MPI.COMM_WORLD)
            else
                MPI.Recv!(untenbek,convert(Int64,id+xprocs),0,MPI.COMM_WORLD)
                x[anzx+2:anzx+2,2:anzy+1]=untenbek
                MPI.Recv!(obenbek,convert(Int64,id-xprocs),0,MPI.COMM_WORLD)
                x[1:1,2:anzy+1]=obenbek
                MPI.Send(unten,convert(Int64,id+xprocs),0,MPI.COMM_WORLD)
                MPI.Send(oben,convert(Int64,id-xprocs),0,MPI.COMM_WORLD)
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
   
   
    anzahl_iterationen=20
    #es wird herausgefynden wie viele prozessore und welche id diese haben
    id = MPI.Comm_rank(MPI.COMM_WORLD)
    num_procs = MPI.Comm_size(MPI.COMM_WORLD)
    
    globx=1
    globy=1
    anzx=0
    anzy=0
    
    #die prozessorenanzahl darf nicht null sein und muss eine Quadratzahl oder 2er Potenz sein(mit bitweise ueberpruefung)
    if num_procs==0 || num_procs&(num_procs-1)!=0 
        if convert(Int64, floor(sqrt(num_procs)))^2!=num_procs 
            return -1
        end
    end
    
    #Die globalen x werte werden bestimmt und die matrixgroesse und das bedingt pro prozessor
    if floor(sqrt(num_procs))^2==num_procs || num_procs==1
        globy=convert(Int64,ceil((n+1)/sqrt(num_procs))*(id%sqrt(num_procs))+1)
        if id%sqrt(num_procs)==sqrt(num_procs)-1
            anzy=convert(Int64,n-globy+1)
        else
            anzy=convert(Int64,ceil((n+1)/sqrt(num_procs)))
        end
        
        globx=convert(Int64,floor(id/sqrt(num_procs))*ceil((n+1)/sqrt(num_procs))+1)
        
        if id>=num_procs-sqrt(num_procs)
            anzx=convert(Int64,n-globx+1)
        else
            anzx=convert(Int64, ceil((n+1)/sqrt(num_procs)))
        end
    else
        xprocs=sqrt(convert(Int64, num_procs/2))
        yprocs=xprocs*2 
        globy=convert(Int64, ceil((n+1)/xprocs)*(id%xprocs)+1)
        anzy=convert(Int64, ceil(n/xprocs))
        if id%xprocs==3
            anzy=anzy-1
        end
        globx=convert(Int64, floor(id/yprocs)*ceil((n+1)/yprocs)+1)
        anzx=convert(Int64, ceil(n/yprocs))
        if floor(id/yprocs)==yprocs-1
            anzx=anzx-1
        end
    end
    
    #Funktionen und A und x werden erzeugt
    f(x,y)=sin(k*pi*x)*sin(k*pi*y)
    f2(x,y)=-2*(k*k*pi*pi)*sin(k*pi*x)*sin(k*pi*y)  
    A=[((n+1)*(n+1)), ((n+1)*(n+1)), -4*((n+1)*(n+1)), ((n+1)*(n+1)), ((n+1)*(n+1))]
    x=fourier(f,n,k,globx,globy,anzx,anzy)
    b=zeros(anzx,anzy)

    #Gauss Seidel
    #x=redblackgs(id,num_procs,x,A,b,anzahl_iterationen)
    x=redblackgs_mit_residuum(id,num_procs,x,A,b,0.1)
    
    #Maximum bestimmen
    max=maximum(abs(x))
    y=MPI.Reduce(max, MPI.MAX, 0, MPI.COMM_WORLD)
    return y
end

MPI.Init()
id = MPI.Comm_rank(MPI.COMM_WORLD)
num_procs = MPI.Comm_size(MPI.COMM_WORLD)
n=127
k=convert(Int64, floor(n/10*4))
if(id==0)
   str= "Daten/red_black_gauss_seidel_$(num_procs).txt"
   datei=open(str,"w")
end
while n<5000
    if(id==0)
        println(n)
        t1=time()
        error=main(n,k)
        t2=time()-t1
        println(n, " ", t2)
        write(datei,string(n),";",string(k),";"string(error),";",string(t2),";","0.1", ";\n")
    else
        error=main(n,k)
    end
    n=2*(n+1)-1
    k=convert(Int64, floor(n/10*4))
end
if(id==0)
    close(datei)
end
MPI.Finalize()
