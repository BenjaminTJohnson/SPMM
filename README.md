# SPMM
Single Particle Melting Model for ice-phase precipitation hydrometeor melting simulations.
Requires CMake 3.20 or newer and a Fortran compiler (tested with gfortran, ifort and ifx).

## Build and run

```
cmake -S . -B build            # add -DCMAKE_BUILD_TYPE=Debug for bounds checking
cmake --build build
mkdir example && cd example
ln -s ../shape_files/shape.dat .
../build/src/spmm
```

To choose a compiler, set `FC` when configuring, for example `FC=ifx cmake -S . -B build`.
`-DSPMM_NATIVE=ON` tunes for the build machine's CPU.

`spmm` (`src/domain_melt5.f90`) reads `shape.dat` and writes one frame per target melt
fraction, named `domain5_fNNNNNN.out`, where `NNNNNN` is the melted mass fraction in
thousandths of a percent. Each line is

```
index  x  y  z  type  body
```

`type` is 1 for ice and -1 for liquid. `body` numbers the connected bodies (26 neighbours,
ice and liquid together); body 1 is the largest, so anything higher is detached. Dipoles are
removed as ice melts so that mass is conserved (ice 916.7, water 999.84 kg m^-3).

## Viewing a run

`tools/melt_viewer.html` is a self-contained 3-D viewer. Open it in a browser, then use
"Open folder" or "Open frames" to load the frame files of any run, including one that is
still going. It needs an internet connection only for fonts and one small library.

`scripts/pack_melt_frames.py OUT.json FRAME...` packs a set of frames into one small JSON
file that the viewer also opens (needs numpy and scipy).

## Surface-energy prototype (`spmm_bond`)

`src/spmm_bond.f90` is an experimental replacement for the liquid movement rule. Liquid
moves to lower a pair energy with separate ice and water affinities, using a noisy phase
followed by a quench. It is not yet a replacement for `spmm`:

- it does not remove dipoles on melting, so it does not conserve mass;
- the fully melted drop ends about 6 cells from the starting centroid, for reasons not yet understood;
- it writes five columns (no body column).

It reads `shape.dat` and an optional namelist file `spmm.nml`. Settings that gave a nearly
spherical final drop (axis ratios 0.988 : 0.994 : 1) on `shape_files/shape.dat`:

```
&spmm
 do_melt=.true., prefix="agg", a_i=1.0,
 temp=4.0, sweeps_hot=10, max_sweeps=40,
 temp_final=8.0, final_hot=10000
/
```

`temp` is the noise level while melting and `temp_final` the level for the final rounding;
levels of 4 or less leave the drop near 0.93 : 0.95 : 1. `a_i` is the ice affinity relative
to water (`a_w`): below 1 the water beads, at 1.0 to 1.2 it coats the ice. `pin_coating`
holds liquid that touches ice in place. Test shapes and analysis scripts are in
`scripts/proto/`.
