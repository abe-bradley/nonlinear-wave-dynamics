# Nonlinear wave dynamics

This is my MATLAB toolkit for simulating nonlinear optical waves and studying them through the lens of statistical mechanics, which we call *optical thermodynamics*.

When you send light into a multimode nonlinear system at moderate power levels, more often than you might expect, the power spreads out among the modes into something that looks thermal. Most of what's here is about watching that happen and checking it against theory.

What you'll find:
- **A lattice builder**: an interactive editor for designing coupled-waveguide lattices, from a single chain to composite lattices with several sublattices. It's my favorite part!
- **Discrete NLSE (DNLSE) simulations** of those lattices, with mode statistics compared against Rayleigh–Jeans and exact microcanonical predictions.
- **Beam propagation (BPM)** codes for the continuous NLSE, in space (2D) and in time (1+1D).
- **Visualizers** for mode dynamics, conservation laws, statistics, and the geometry of phase space.

More is on the way: a multimode NLSE solver and Lugiato–Lefever codes (see [Coming soon](#coming-soon)). This is active research code, so things will keep changing.


## Getting started

1. Clone or download this repository.
2. In MATLAB, go to the repository folder and run `init_dir`. It adds the repository to your MATLAB path and moves you into its top folder. Heads up: it also clears your workspace.
3. Open any script in `MATLAB/Operators` and run it.

The scripts in `MATLAB/Operators` are starting points, in the sense that each one has its parameters near the top, with lots of alternatives (other lattices, other launch conditions, other structures) sitting right there commented out. Changing parameters and uncommenting things is a great way to get familiar. The defaults take anywhere from a few seconds to about a minute to run.

**Requirements.** So far, just MATLAB, with a desktop session for the interactive tools. I develop in R2026b, and I have tested most of the repo on version dating back to R2024a. The Symbolic Math Toolbox is optional: the microcanonical code only reaches for it when double precision can't resolve a cancellation, and it tells you if that happens.

## The lattice builder

```matlab
[H, B] = lattice_config(true);
```

This opens a full graphical editor for designing a lattice, then hands you its Hamiltonian `H`, ready to drop into a simulation.

**Designing a sublattice**

- **Start from a preset or from scratch.** Square, triangular, honeycomb, Lieb, ring, and the two-dimensional $J_x$ lattice (couplings proportional to $\sqrt{n(N-n)}$ along each axis, which gives an evenly spaced spectrum). Or pick Custom and place every site yourself.
- **Edit directly on the canvas.** Double-click to add a site, drag to move, right-click to delete, Shift-click or box-select to grab several at once, and Ctrl+Z to undo. Scroll to zoom and use the hand button to pan.
- **Couplings come from geometry.** Pick a distance law, a decay parameter $\gamma$, and a range cutoff, and every pair of sites gets its coupling automatically. With $q = d/d_{\min}$ (distance over the shortest site spacing), the laws are
  ```math
  \kappa_{\text{power}} = q^{-\gamma}, \qquad
  \kappa_{\text{exp}} = e^{-\gamma (q-1)}, \qquad
  \kappa_{\text{fiber}} = q^{-1/2}\, e^{-\gamma (q-1)}.
  ```
- **Override any coupling by hand.** Click an edge and set its magnitude and phase. Each direction is independent, so couplings don't have to be reciprocal. "Restore geometric couplings" puts everything back.
- **Detune individual sites.** Give nodes on-site detunings, with disorder sitting on top.
- **Disorder.** Several common forms of disorder are supported, namely Anderson, coupling, or Aubry–André (quasiperiodicity). Disorder samples stay attached to their sites while you keep editing.

**Assembling a composite lattice**

Add as many sublattices as you like, each with its own internal design, then drag them around the editing plane or type in positions and rotations. Couplings *between* sublattices follow their own distance law, decay, and cutoff, while each sublattice's internal couplings stay exactly as you designed them.

The output `H` is ordered sublattice by sublattice, and `B{k}` lists the sites of sublattice `k`, so `H(B{k},B{k})` is its block. The simulations use `B` to compute statistics in each sublattice's own modes as well as for the full lattice.

**Saving, loading, and viewing modes**

**Save** writes the whole editable design to a `.mat` file (by default into the example-lattice folder), and **Load** brings it back for more editing. **Build** returns `H` and `B`. In your own scripts you can load a saved lattice without opening the editor at all:

```matlab
lattices = fullfile(codes_root(), 'MATLAB', 'Data', ...
    'System configurations', 'Lattice configurations');
[H, B] = lattice_config(true, fullfile(lattices, 'Jx_4x4'));
```

With the first argument set to `true`, a mode viewer opens afterward, showing the spectrum and every mode's intensity and phase drawn on the lattice itself, for the full lattice and for each sublattice.

**Example lattices** (in `MATLAB/Data/System configurations/Lattice configurations`):

| File | Sites | Description |
| --- | --- | --- |
| `4_line`, `5_line` | 4, 5 | Short chains |
| `14_ring` | 14 | Ring |
| `Jx_3x2` to `Jx_6x6` | 6 to 36 | Two-dimensional J<sub>x</sub> lattices |
| `double_test` | 9 + 3 | Composite lattice with two sublattices |
| `triple_test` | 9 + 7 + 8 | Composite lattice with three sublattices |

## Simulating lattices: the DNLSE

The lattice simulations solve the discrete nonlinear Schrödinger equation

```math
i\,\frac{da_n}{dz} = \sum_m H_{nm}\, a_m + |a_n|^2 a_n ,
```

where $a_n$ is the field in waveguide $n$ and $H$ comes from the lattice builder. For a Hermitian $H$, two quantities are conserved: the total power $P = \sum_n |a_n|^2$ and the Hamiltonian $\sum_{nm} a_n^{\ast} H_{nm} a_m + \tfrac{1}{2}\sum_n |a_n|^4$.

`MATLAB/Operators/Propagation/DNLSE/kerr_dnlse.m` runs the whole sequence:

1. **Load a lattice** (`Jx_4x4` by default) and open its mode viewer.
2. **Launch an ensemble.** Spread some power across a few chosen modes, then add amplitude and phase noise and randomize the phases, so every member of the ensemble is a different realization.
3. **Propagate** each member with `ode45` at tight tolerances.
4. **Project onto the modes** to get the occupations $J_m = |c_m|^2$, for the full lattice and for each sublattice.
5. **Check conservation laws.** A viewer shows the power and the linear, nonlinear, and total energy. It also splits the nonlinear energy into a uniform-phase part and a four-wave-mixing part.
6. **Collect statistics** over the late portion of propagation and across all members, and compare them with theory.

Two more viewers open along the way. The **modal dynamics** viewer shows mode evolution, entropy evolution, three-mode trajectories, complex-plane portraits of each mode, and mutual coherence. The **thermodynamics** viewer compares simulation against theory: mean occupations, reciprocal occupations, marginal distributions, Pearson correlations, and two-mode joint distributions.

### The theory side

Both predictions use the power $P$ and the linear energy $U = \sum_m \varepsilon_m J_m$, where the $\varepsilon_m$ are the eigenvalues of $H$. This is the natural setting when the nonlinearity is weak compared to the coupling.

**Rayleigh–Jeans.** This is the standard optical thermodynamics result,

```math
\langle J_m \rangle = \frac{T}{\varepsilon_m - \mu},
```

with the temperature $T$ and chemical potential $\mu$ fixed by $P$ and $U$.

**Exact microcanonical.** For a finite number of modes $M$, I count every arrangement of occupations that has *exactly* power $P$ and energy $U$. In available phase-space volume is then

```math
\Omega(U,P) = \int_{J_m \ge 0} \delta\Big(P - \sum_m J_m\Big)\, \delta\Big(U - \sum_m \varepsilon_m J_m\Big)\, d^M J ,
```

which the code evaluates in closed form, including degenerate spectra. The expressions can get computationally complex rather rapidly as the system scales in size, so the numerical implementations are designed carefully enough to avoid catastrophic cancellations. Everything else follows from it: the temperature $T = (\partial \ln\Omega / \partial U)^{-1}$, the chemical potential, the mean occupations, the full marginal distributions, and the correlations between modes. With only a handful of modes the average modal content differs noticeably from Rayleigh–Jeans, and the two generally agree as $M$ grows.

These functions also work on their own:

```matlab
[T, mu, J] = calculate_MC(U, P, epsilon);   % exact microcanonical
[T, mu, J] = calculate_RJ(U, P, epsilon);   % Rayleigh–Jeans
```

## Seeing the constraints

Fixing $P$ and $U$ confines the occupations to a convex polytope: all nonnegative occupations that add up to $P$, sliced by the energy plane $\sum_m \varepsilon_m J_m = U$. Two scripts in `MATLAB/Operators/Misc plots` let you look at it directly.

- `constrained_dynamics.m` runs the DNLSE pipeline on `Jx_3x3`, then animates the trajectories of $(J_1, J_2, J_3)$ moving through the allowed region, side by side with the relative phases. You can play, pause, scrub, rotate, and switch to the ensemble average.
- `manifold_visualization.m` is an interactive explorer of that polytope. Choose which three occupations to plot and use sliders for the rest. The launch point and the microcanonical equilibrium are marked.

## Beam propagation (BPM)

Both scripts in `MATLAB/Operators/Propagation/BPM` use the split-step Fourier method.

- **`bpm_xz.m`: beams in space.** A paraxial beam in a 2D refractive-index landscape $n(x,z)$ with Kerr nonlinearity. A built-in finite-difference mode solver finds the guided modes, so you can launch them directly. The default is a single waveguide with its first two modes launched. Commented out and ready to go: a single soliton, two interacting solitons with a relative phase, a waveguide with a Gaussian launch, a directional coupler, and manual Gaussian or sech launches. The index is stored at every $z$ step, so structures that change along $z$ just mean filling in that array.
- **`bpm_tz.m`: light in time.** The classic fiber NLSE, with dispersion from a fused-silica Sellmeier model, loss, and self-phase modulation. It handles pulses like any temporal split-step code. The configuration I pushed, though, is a finite-linewidth CW laser: an ensemble of fields with random-walk phase noise, whose averaged spectra show how the line broadens. For classic pulse propagation, replace the initial field with your pulse and use a single ensemble member.

## Repository layout

```
init_dir.m        adds the repository to the MATLAB path
MATLAB/
  Operators/      the scripts you run
  Functions/      the building blocks they call
  Data/           example lattices
```

Inside `MATLAB/Functions`:

| Folder | Contents |
| --- | --- |
| `System configuration` | `lattice_config`, the lattice builder |
| `Propagators` | `dnlse`, the DNLSE right-hand side |
| `Mode solvers` | `eigk` (modes plus Kerr overlap coefficients), `laplacian_1d` |
| `Conserved quantities` | `kerr_conservation` |
| `Statistics` | `modal_statistics`: moments, correlations, histograms |
| `Microcanonical gas` | `volume`, `mixed_moments`, `calculate_MC`, `calculate_RJ`, `marginals`, `pearson_correlations` |
| `Plotting` | the viewers, a progress bar, and `ot_colormap` |
| `Materials` | `sellmeier_silica` |
| `Basics` | `codes_root` (repository location), `noiser`, `degeneracy`, and small helpers |

`ot_colormap` gives every plot one interface to MATLAB's built-in colormaps and about 180 bundled scientific palettes.

## Coming soon

These are on their way. No dates, but they're coming.

- **Multimode NLSE (MMNLSE).** A coupled-mode solver for multimode fibers and waveguides. It propagates the complex amplitudes of the guided modes, coupled to each other through four-wave mixing and random linear mode coupling. The propagation engine is basically done, and it's blazing fast. The tools for setting up the excitation are what I'm still building. Mode solving with COMSOL, including the COMSOL models, will arrive with it. It plugs straight into the same statistics and thermodynamics viewers described above.
- **Lugiato–Lefever equation (LLE).** Codes for driven, lossy Kerr resonators. Turns out they thermalize too!

## License

MIT. See [`LICENSE`](LICENSE).
