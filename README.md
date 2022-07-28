# SPMM
Single Particle Melting Model for ice-phase precipitation hydrometeor melting simulations 
Requires any reasonably modern fortran compiler.

Instructions:
```
mkdir build
cd build
cmake ..
make
cd ..
mkdir example
cd example
ln -s ../shape_files/shape.dat .
ln -s ../build/src/spmm .
./spmm
```
