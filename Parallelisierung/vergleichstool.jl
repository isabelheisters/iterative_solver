#NIEMALS DIE REFERENZ VERAENDERN
referenz = open("Daten/matrixmultiplikation_4_threads.txt", "r")
#akt ist die Datei mit der man vergleichen soll
akt = open("Daten/matrixmultiplikation_8_threads.txt", "r")

# Datei sieht folgendermassen aus:
# N, k,fehler
vergleich=readlines(referenz)
i=1
j=1
aktuell=readlines(akt)
vline=split(aktuell[i],";")
aktline=split(aktuell[j], ";")
close(referenz)
close(akt)

koeffizient=0.0001

while(parse(Int,aktline[1])<parse(Int,vline[1]))
    j=j+1
    aktline=split(aktuell[j], ";")
end
while(parse(Int,aktline[1])>parse(Int,vline[1]))
    i=i+1
    vline=split(aktuell[i],";")
end

while(i<=size(vergleich)[1]&&j<=size(aktuell)[1])
    if(parse(Int, aktline[1])==parse(Int,vline[1]))
        print("Vergleichswert: ", aktline[1])
        if(parse(Float64,vline[3])-parse(Float64,aktline[3])<= koeffizient)
            println(" --> Stimmt Fehler: ",string(parse(Float64,vline[3])-parse(Float64,aktline[3])),"\n")
            i=i+1
            j=j+1
            if(i<=size(vergleich)[1]&&j<=size(aktuell)[1])
                vline=split(aktuell[i],";")
                aktline=split(aktuell[j], ";")
            end
        else
            print(" --> Stimmt nicht... Fehler: ",parse(Float64,vline[3])-parse(Float64,aktline[3]),"\n")
            i=Nan
        end
    else
        println("Du hast anscheinend nicht mit den richtigen Ns gearbeitet\n")
    end
end
