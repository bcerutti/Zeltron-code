# Remove all .mod and .o files
rm *.mod *.o 
rm zeltron.exe zeltron.out zeltron.err rankLog_*.txt

# Remove all data
rm -r ./data/*
rm ./data_restore/*

if [ ! -d data ]
then
  mkdir data
fi

# Create folders where data will be dumped
mkdir data/orbits
mkdir data/fields
mkdir data/densities
mkdir data/currents
mkdir data/particles
mkdir data_restore

# Just make the executable; don't run it.
make IO=txt

# Run the code
mpirun -np 1 zeltron.exe
