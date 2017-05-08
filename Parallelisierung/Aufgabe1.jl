function zerlegung_mit_pivot(A)
    P=ones(size(A)[1])
    for k in 1:size(P)[1]
        P[k]=k
    end
    for i in 1:size(A)[1]-1
        maxzeile=i
        for m in i+1:size(A)[1]
            if(A[maxzeile,i]<A[maxzeile,m])
                maxzeile=m
            end
        end
        hilf=A[maxzeile:maxzeile,1:end]
        A[maxzeile:maxzeile,1:end]=A[i:i,1:end]
        A[i:i,1:end]=hilf
        hilf=P[maxzeile]
        P[maxzeile]=P[i]
        P[i]=hilf    
        for k in i:size(A)[2]
            faktor=A[k,i]/A[i,i]
            for j in i:size(A)[2]
                if i<j
                    if k!=i
                        A[k,j]=A[k,j]-faktor*A[i,j]
                    end
                elseif k!=j&&k!=i
                    A[k,j]=A[k,j]/A[i,j]
                end
            end
        end
    end
    return A,P
end

function permutation(p,x)
    y=zeros(size(p)[1])
    for i in 1:size(p)[1]
        hilf=convert(Int64,p[i])
        y[i]=x[hilf]
    end
    return y
end

function vorwaerts(lu,b)
    x=zeros(size(A)[1])
    for i in range(1,1,size(lu)[1])
        x[i] = b[i]
        if i!=size(lu)[1]
            for j in 1:i
                b[i+1]=b[i+1]-lu[i+1,j]*x[j]
            end
        end 
    end
    return x
end

function rueckwaerts(lu,b)
    x=zeros(size(A)[1])
    for i in range(size(A)[1],-1,size(A)[1])
        x[i] = b[i]/A[i,i]
        if i!=1
            for j in i:size(A)[1]
                b[i-1]=b[i-1]-A[i-1,j]*x[j]
            end
        end     
    end
    return x
end

function routine(A,b)
    lu,p=zerlegung_mit_pivot(A)
    b=permutation(p,b)
    b=vorwaerts(lu,b)
    b=rueckwaerts(lu,b)
    return b
end

A=[0.0 0 0 1; 2 1 2 0;4 4 0 0;2 3 1 0]
b=[3,5,4,5]
routine(A,b)
println(A)
A=[0.0 0 0 1; 2 1 2 0;4 4 0 0;2 3 1 0]
b=[4,10,12,11]
routine(A,b)
