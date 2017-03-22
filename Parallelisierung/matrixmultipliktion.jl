using MPI

function main(n::Int64, k::Int64)
    #k=1

    id = MPI.Comm_rank(MPI.COMM_WORLD)
    size = MPI.Comm_size(MPI.COMM_WORLD)
    #n=256
    if(n%size==0)
        lokn=convert(Int64,(n/size))
    elseif(n%size==size-1)
        if(id==size-1)
            lokn=convert(Int64,(n/size)-1)
        else
            lokn=convert(Int64,(n/size))
        end
    else
        MPI.Finalize()
        return 0
    end
    
    f(x)=sin(k*pi*x)
    f2(x)=-((k*pi)^2)*sin(k*pi*x)   
    A=[((n+1)*(n+1)), -2*((n+1)*(n+1)), ((n+1)*(n+1))]
    x=fourier(f,lokn,k,n,id)
    x=fuege_Rand_hinzu(x)
    
    #Hier kommt die Kommunikation hin!
    w=zeros(1)
    if(id%2==0)
        if(id==0)
            MPI.Send([x[lokn+1]],id+1,0,MPI.COMM_WORLD)
            MPI.Recv!(w,id+1,0,MPI.COMM_WORLD)
            x[lokn+2]=w[1]
        elseif(id==size-1)
            MPI.Send([x[2]],id-1,0,MPI.COMM_WORLD)
            MPI.Recv!(w,id-1,0,MPI.COMM_WORLD)
            x[1]=w[1]
        else
            MPI.Send([x[2]],id-1,0,MPI.COMM_WORLD)
            MPI.Recv!(w,id-1,0,MPI.COMM_WORLD)
            x[1]=w[1]
            MPI.Send([x[lokn+1]],id+1,0,MPI.COMM_WORLD)
            MPI.Recv!(w,id+1,0,MPI.COMM_WORLD)
            x[lokn+2]=w[1]
        end
    else
        if(id==size-1)
            MPI.Recv!(w,id-1,0,MPI.COMM_WORLD)
            x[1]=w[1]
            MPI.Send([x[2]],id-1,0,MPI.COMM_WORLD)
        else
            MPI.Recv!(w,id-1,0,MPI.COMM_WORLD)
            x[1]=w[1]
            MPI.Send([x[2]],id-1,0,MPI.COMM_WORLD)
            MPI.Recv!(w,id+1,0,MPI.COMM_WORLD)
            x[lokn+2]=w[1]
            MPI.Send([x[lokn+1]],id+1,0,MPI.COMM_WORLD)
        end
    end

    x=matrixmultiplikation(A,x)

    
    x=norm(fourier(f2,lokn,k,n,id)-x,Inf)

    y=MPI.Reduce(x, MPI.MAX, 0, MPI.COMM_WORLD)
    if(id==0)
        println("Error: ",y)
    end
    return y

end



function fuege_Rand_hinzu(u::Array{Float64,1})
    N=size(u)[1]
    matrix=zeros(N+2)
    matrix[2:N+1] = u
    return matrix
end

function matrixmultiplikation(A::Array{Int64,1},x::Array{Float64,1})
    gr=size(x)[1]
    matrix=zeros(gr-2)
    for i in 2:gr-1
        matrix[i-1]=A[1]*x[i-1]+A[2]*x[i]+A[3]*x[i+1]
    end
    return matrix
end

function fourier(f::Function, N::Int64,k::Int64, Nglob::Int64, akt::Int64)#Nsollte auch nur die lokale Anzahl sein 
#nicht die globale
    matrix=zeros(N)
    l=1
    for i in akt*N+1:(akt+1)*N
        matrix[l]=f(i/(Nglob+1)) 
        l=l+1
    end
    return matrix
end

MPI.Init()
id = MPI.Comm_rank(MPI.COMM_WORLD)
n=40
k=38
if(id==0)
    datei=open("Daten/matrixmultiplikation_8_threads.txt","w")
end
while n<10000
    error=main(n,k)
    if(id==0)
        write(datei,string(n),";",string(k),";"string(error), ";\n")
    end
    n+=8
end
if(id==0)
    close(datei)
end
MPI.Finalize()
