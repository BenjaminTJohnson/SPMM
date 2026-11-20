# SPMM
Single Particle Melting Model for ice-phase precipitation hydrometeor melting simulations.
It takes an ice particle described as dipoles on a cubic lattice and melts it step by step,
writing the particle at a series of melted fractions for use in scattering calculations.

Requires CMake 3.20 or newer and a Fortran compiler (tested with gfortran, ifort and ifx).
The viewer and helper scripts need Python 3 with numpy and scipy.

## Build and run

```
cmake -S . -B build            # add -DCMAKE_BUILD_TYPE=Debug for bounds checking
cmake --build build
mkdir example && cd example
ln -s ../shape_files/shape.dat .
cp ../shape_files/spmm.nml .
../build/src/spmm
```

To choose a compiler, set `FC` when configuring, for example `FC=ifx cmake -S . -B build`.
`-DSPMM_NATIVE=ON` tunes for the build machine's CPU.

The sample takes about 3.5 minutes and writes 26 frames, `agg_f000000.out` to
`agg_f100000.out`. It starts with 159,957 dipoles and ends with 146,657 (13,300 removed for
mass conservation) as a single, nearly spherical drop with axis ratios 0.988 : 0.994 : 1.

`python3 ../scripts/run_packed.py` does the same run and also packs the frames into one
small file for the viewer (see "Viewing a run").

## Input

`spmm` reads `shape.dat` from the run folder: one header line, then one dipole per
line as `x y z type`, with integer lattice coordinates. `type` 2 is liquid; anything else is
ice. Settings are read from an optional namelist file `spmm.nml` in the same folder.

Use `shape_files/spmm.nml` as the starting point. Without it the noise levels default to
zero and the fully melted drop stays visibly faceted (near 0.93 : 0.95 : 1).

## Output

Frames are named `<prefix>_fNNNNNN.out`, where `NNNNNN` is the melted mass fraction in
thousandths of a percent. Each line is

```
index  x  y  z  type  body
```

`type` is 1 for ice and -1 for liquid. `body` numbers the connected bodies (26 neighbours,
ice and liquid together); body 1 is the largest, so anything higher is detached.

A frame is written at 0%, at each target melted fraction, and at 100%. The default targets
are 0.5%, every 1% from 1% to 10%, every 2.5% to 20%, then 25%, 30%, every 10% to 90%, and
95%: 26 frames, 16 of them at or below 20%, where absorption changes fastest. Each melt step
melts just enough ice to reach the next target, so frames land on the targets to within one
dipole, but never more than a cap: 0.5% of the particle per step below 20% melted, 2% above.
Which eligible ice melts in a step is chosen at random from a fixed seed.

To change this, add a `&melt` group to `spmm.nml`:

```
&melt
 targets = 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.05, 0.10, 0.20, 0.50,
 dfmax = 0.005, dfmax_late = 0.02, f_late = 0.20, seed = 12345
/
```

`targets` are fractions between 0 and 1 in ascending order (up to 100 of them); `dfmax` and
`dfmax_late` are the caps below and from `f_late`.

## How it works

Ice melts where it is exposed: an ice dipole with few ice neighbours is eligible, and the
threshold rises when nothing qualifies. Liquid then moves to lower a pair energy with
separate ice and water affinities, in a noisy phase followed by a quench. The `&spmm` group
of `spmm.nml` controls this:

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

Four things are done by default, each with a switch:

- **Mass conservation.** Ice (916.7 kg m^-3) melts to a smaller volume of water (999.84), so
  about one melting dipole in twelve is removed instead of becoming water. The melted
  fraction is the melted mass fraction, counting removed dipoles as melted.
  `conserve_mass=.false.` keeps every dipole.
- **Reattaching detached pieces.** Each water body only rearranges itself, so pieces can
  separate: a water blob, or an ice fragment whose connecting neck has melted and drained.
  Left alone they never rejoin (on a 20-flake aggregate the fully melted particle ended as
  104 separate drops). After every melt step each body that has lost contact with the rest
  is shifted rigidly, without changing shape, to the nearest position where it touches
  another body; the largest body never moves. Water that touches other water then merges
  with it in the normal relaxation. `reattach=.false.` turns this off, for diagnosis only.
- **Centre-of-mass frame.** Ice is held on the lattice while water relocates toward the bulk
  of the particle, so the centre of mass wanders in lattice coordinates (by 5 to 9 cells
  over a full melt), which a free particle's cannot. Each frame is written shifted by whole
  cells so that its centre of mass stays within half a cell of where it started; the
  simulation itself is not changed. `recentre=.false.` writes lattice coordinates.

- **Minimum water thickness.** Discrete-dipole scattering codes misrepresent water that is
  only one dipole thick (absorption too low by a third or more at Ku and Ka band), and
  melting starts as scattered single water dipoles and one-dipole films. Before each frame
  is written the water is gathered until every water dipole lies in a 2 x 2 x 2 cube of
  water: thin water next to a nearly complete cube moves in to complete it, nearest first,
  travelling along the particle surface. A few dipoles (about ten at most on the test below)
  can be left over where no cube can be completed, and the first frames, with only a handful
  of cubes, move water furthest. On a 12,708-dipole test particle 95 to 100% of the water
  lies in such cubes in every frame, against 0% at 2% melted and 48% at 10% melted without
  the rule. `water_thick` sets the cube size (default 2); `water_thick=1` turns it off.

Test shapes and analysis scripts are in `scripts/proto/`;
`python3 scripts/proto/axes.py FRAME...` prints axis ratios, body counts and the centroid.

## Viewing a run

`tools/melt_viewer.html` is a self-contained 3-D viewer. Open it in a browser, then use
"Open folder" or "Open frames" to load the frame files of any run, including one that is
still going. It needs an internet connection only for fonts and one small library.

`scripts/pack_melt_frames.py OUT.json FRAME...` packs a set of frames into one small JSON
file that the viewer also opens (a few megabytes instead of hundreds).

`scripts/run_packed.py` runs `spmm` in a run folder and packs the frames into
`packed.json` there when it finishes:

```
cd my_run                               # holds shape.dat and spmm.nml
python3 ../scripts/run_packed.py        # add --live to repack as each frame is written
```

With `--live` the JSON is rebuilt whenever a new frame is complete, so it can be opened and
reopened in the viewer during the run. `--exe PATH` runs a different program, `--dir` names
the run folder, and `--clean` deletes frames left there by an earlier run (without it the
script stops if it finds any).

## DDSCAT targets

`scripts/frames_to_ddscat.py` turns frames into DDSCAT target files in the `FROM_FILE`
format (needs numpy):

```
python3 scripts/frames_to_ddscat.py --spacing 20 my_run/agg_f*.out
```

This writes `my_run/ddscat/<frame>/shape.dat`, one folder per frame, and a table
`my_run/ddscat/frames.txt` with each frame's melted fraction, dipole counts and effective
radius. Composition 1 is ice and composition 2 is liquid water in every file, so
`ddscat.par` needs `NCOMP = 2` with the ice refractive-index file first and the water file
second, including for frames that are all ice or all water. DDSCAT takes the dipole spacing
from the effective radius, `aeff = (3 N / 4 pi)^(1/3) d`: the table gives `aeff/d`, and with
`--spacing` (the dipole spacing in micrometres) also `aeff` in micrometres. `--out` names
the output folder and `--main-body` keeps only the largest connected body.

## Deprecated: `spmm_legacy`

The earlier program (`src/domain_melt5.f90`) is deprecated and is not built by default.
Configure with `-DSPMM_BUILD_LEGACY=ON` to build it as `spmm_legacy`. It uses the same
melting rule, input and frame format, with frames named `domain5_fNNNNNN.out`, and reads
only the `&melt` group of `spmm.nml`. Its liquid moves one cell at a time toward the
particle's centre of mass, and only when it is not touching ice. Compared with `spmm`:

- the fully melted drop is less round (0.835 : 0.927 : 1 on the sample shape);
- on open aggregates pieces detach in mid-melt and are not brought back;
- frames are in lattice coordinates, with a drift of 1 to 3 cells.

It conserves mass in the same way and ends with the same number of dipoles.
