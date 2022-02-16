#/bin/sh
# Refresh parameters

File='settings/parameters.txt'

rm -rf Data
mkdir -p Data

for ii in 1 2 3 4 5
do 
    for jj in 1 2 3 4 5
    do 
        Rmin=$(sed -n 3p settings/parameters.txt | tr -d -c 0-9.-)
        sed -i "26s/.*/$ii          #nDeltaT/" $File
        sed -i "27s/.*/$jj          #nEpsilonT/" $File
        make
        make plot
        cp results/data.txt Data/d_${ii}e_${jj}.txt
        rm -rf results/graph
    done
done
