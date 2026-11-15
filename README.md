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

`shape.dat` has one header line, then one dipole per line as `x y z type`.

## Output frames and melt steps

Both programs write a frame at 0%, at each target melted fraction, and at 100%. The default
targets are 0.5%, every 1% from 1% to 10%, every 2.5% to 20%, then 25%, 30%, every 10% to
90%, and 95%: 26 frames, 16 of them at or below 20%, where absorption changes fastest. Each
melt step melts just enough ice to reach the next target, so frames land on the targets to
within one dipole, but never more than a cap: 0.5% of the particle per step below 20%
melted, 2% above. Which eligible ice melts in a step is chosen at random from a fixed seed.

To change this, put a `&melt` group in a file `spmm.nml` in the run folder (read by both
programs):

```
&melt
 targets = 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.05, 0.10, 0.20, 0.50,
 dfmax = 0.005, dfmax_late = 0.02, f_late = 0.20, seed = 12345
/
```

`targets` are fractions between 0 and 1 in ascending order (up to 100 of them); `dfmax` and
`dfmax_late` are the caps below and from `f_late`.

## Viewing a run

`tools/melt_viewer.html` is a self-contained 3-D viewer. Open it in a browser, then use
"Open folder" or "Open frames" to load the frame files of any run, including one that is
still going. It needs an internet connection only for fonts and one small library.

`scripts/pack_melt_frames.py OUT.json FRAME...` packs a set of frames into one small JSON
file that the viewer also opens (needs numpy and scipy).

## Surface-energy prototype (`spmm_bond`)

`src/spmm_bond.f90` replaces the liquid movement rule. Liquid moves to lower a pair energy
with separate ice and water affinities, using a noisy phase followed by a quench. It writes
frames in the same six-column format as `spmm`, named `<prefix>_fNNNNNN.out`, and conserves
mass in the same way (`conserve_mass=.false.` turns the removal of dipoles off). One thing is
not yet understood: the fully melted drop ends about 5 cells from the starting centroid.

It reads `shape.dat` and an optional namelist file `spmm.nml`. Settings that give a nearly
spherical final drop (axis ratios 0.982 : 0.994 : 1, 146,657 dipoles, about 3.5 minutes) on
`shape_files/shape.dat`:

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
holds liquid that touches ice in place.

Each water body only rearranges itself, so pieces can separate: a water blob, or an ice
fragment whose connecting neck has melted and drained. Left alone they never rejoin (on a
20-flake aggregate the fully melted particle ended as 104 separate drops). So after every
melt step each body that has lost contact with the rest is shifted rigidly, without changing
shape, to the nearest position where it touches another body; the largest body never moves.
Water that touches other water then merges with it in the normal relaxation, and the run
prints the number of shifts made. `reattach=.false.` turns this off, for diagnosis only.

Test shapes and analysis scripts are in
`scripts/proto/`.
