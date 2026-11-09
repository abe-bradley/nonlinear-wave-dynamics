# Nonlinear wave dynamics

MATLAB research examples for Kerr nonlinear lattice dynamics, optical thermodynamics, constrained occupations, and beam propagation. This is actively developed research code; numerical choices and interfaces may change.

## Getting started

1. Download or clone the repository into any folder.
2. Open MATLAB and make the folder containing `init_dir.m` your current folder.
3. Run `init_dir`.
4. Open and run one of the examples below. Run initialization first: the example scripts clear workspace variables but keep the MATLAB search path.

The root folder can have any name. Initialization adds that folder and its descendants, excluding `Packages`, `Recycling`, and `.temp`. No parent folder or external package directory is needed for these examples.

## Examples

| Script | Purpose |
| --- | --- |
| `MATLAB/Operators/Propagation/DNLSE/kerr_dnlse.m` | Propagate an ensemble on a nonlinear lattice and compare modal statistics with thermodynamic predictions. Opens the lattice configuration editor. |
| `MATLAB/Operators/Misc plots/constrained_dynamics.m` | Explore constrained modal dynamics and equilibrium statistics. |
| `MATLAB/Operators/Misc plots/manifold_visualization.m` | Interactively explore nonnegative modal occupations at fixed power and energy. |
| `MATLAB/Operators/Propagation/BPM/bpm_xz.m` | Beam propagation in space. |
| `MATLAB/Operators/Propagation/BPM/bpm_tz.m` | Pulse propagation in time. |

All ten example lattices are in `MATLAB/Data/System configurations/Lattice configurations`. Supporting functions retain their original directory organization under `MATLAB/Functions`.

The scripts contain their parameters near the top. Some defaults involve large grids or ensembles and may take substantial time or memory. Reduce propagation grid sizes and ensemble counts when first exploring them; keep dependent parameters consistent (for example, `N_end <= Nz`).

## Requirements and validation

Development and the colormap tests have used MATLAB R2024b. Other releases have not been systematically checked. Interactive examples need MATLAB graphics support.

Some microcanonical calculations can use Symbolic Math Toolbox (`sym` and `vpa`) to recover from numerical cancellation. These branches report an error when that toolbox is unavailable. Available native MATLAB colormaps also depend on the MATLAB release.

Run the colormap tests after initialization:

```matlab
results = runtests(fullfile(codes_root(), 'MATLAB', 'Functions', ...
    'Plotting', 'Colormaps', 'tests'));
assert(all([results.Passed]))
```

These tests check palette selection and sampling. They do not establish the physical accuracy of the simulation examples.

## Colormaps and third-party material

`ot_colormap` selects supported native MATLAB maps or bundled external palettes. For example:

```matlab
colormap(ot_colormap('parula', 256))
colormap(ot_colormap('Set1', 3, 'Sampling', 'colorbrewer'))
```

See `MATLAB/Functions/Plotting/Colormaps/COLORMAPS.txt` for behavior and source provenance. External palette permissions and required credits are in that folder's `THIRD_PARTY_NOTICES.txt` and `LICENSES` directory. Retain these when redistributing the bundled palettes. Native MATLAB maps are requested from the user's MATLAB installation; their RGB tables are not bundled.

## License

The project's original code is available under the MIT License; see `LICENSE`. Bundled third-party material retains its own licenses and attribution requirements, documented in the colormap folder above.
