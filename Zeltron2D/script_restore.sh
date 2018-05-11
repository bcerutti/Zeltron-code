# Remove all .mod and .o files
rm *.mod *.o 

# Compile the code
make IO=txt

# Run the code
mpirun -np 1 zeltron.exe
