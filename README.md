# Wingbeat-Call Asynchrony Simulation and Theory

This repository contains MATLAB code used to simulate and analyse temporal feasibility constraints on wingbeat-call coordination in echolocating bats.

This code supports the following manuscripts:

1. [*Biosonar Responsivity Sets the Stage for the Terminal Buzz* ](https://doi.org/10.1101/2025.06.16.659925)
2. [*Temporal Feasibility Constraints on Wingbeat-Call Synchrony in Actively Echolocating Bats*](https://doi.org/10.1101/2025.06.18.660328)
3. [*Swarm Cohesion in Bats Emerges from Stable Temporal Loops*](https://www.biorxiv.org/content/10.1101/2025.07.05.663265v1)

---

## Repository Structure

### `init.m`
- Adds all subfolders under `repo/` to the MATLAB path.

### `mat/` (main scripts)
- `wingbeat_sim_loop.m`
	- Runs three deterministic scenario simulations:
		- fixed wingbeat frequency, fixed excursion
		- dynamic wingbeat frequency, fixed excursion
		- dynamic wingbeat frequency, dynamic excursion
	- Exports synchrony and waveform figures (PDF).
- `sim_rand_var_run.m`
	- Runs Monte Carlo simulations with randomised parameter sampling.
	- Produces a `results` table and saves `simulation_results.mat`.
- `analyse_sim_results.m`
	- Loads and analyses simulation output.
	- Computes descriptive summaries, effect sizes, and circular statistics.
	- Exports summary CSV files and analysis figures.

### `src/` (helper functions)
- Core simulation and plotting helpers called by scripts in `mat/`.
- Includes `simulateEcholocationWings.m`, `generateVirtualBatCall.m`, and formatting helpers (`formatLatex*.m`).

### `data/`
- Input and derived data files used by the analysis pipeline.
- Includes `simulation_results.mat` and exported summary CSVs.

### `results/`
- Output artefacts from simulation and analysis runs.

### `doc/`
- Supplementary documentation and LaTeX assets.

### Other files
- `citation.cff`: citation metadata.
- `LICENSE.md`: repository licence.

---

## Requirements

- MATLAB R2020b or newer.

---

## Quick Start

From MATLAB, with the working directory set to `repo/`:

```matlab
init
```

Then run scripts in this order:

1. Deterministic scenario figures:

```matlab
run('mat/wingbeat_sim_loop.m')
```

2. Monte Carlo simulations:

```matlab
run('mat/sim_rand_var_run.m')
```

3. Analysis of Monte Carlo outputs:

```matlab
run('mat/analyse_sim_results.m')
```

Notes:
- `sim_rand_var_run.m` writes `simulation_results.mat` in the current working directory.
- `analyse_sim_results.m` expects `../data/simulation_results.mat` by default when `results` is not already in workspace. Update the path in the script if your output location differs.
- Some scripts currently contain absolute output paths for figure export; adjust those paths for your local setup.

---

## Reproducibility Notes

- Monte Carlo runs use random sampling (`rng('shuffle')` in `sim_rand_var_run.m`), so outputs differ between runs unless you set a fixed RNG seed.
- Analysis outputs are exported as CSV tables for direct manuscript reporting and cross-checking.

---

## License

This project is licensed under the [Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)](https://creativecommons.org/licenses/by-nc/4.0/).

You are free to share and adapt the material for non-commercial use with appropriate credit.

## Citation

If you use this codebase in part or in full, please cite the corresponding manuscripts or contact the author for further information.

## Other Publications

Check my [ResearchGate](https://www.researchgate.net/profile/Ravi-Umadi-3) or [ORCID](https://orcid.org/0000-0003-3867-1769) for a full list of my research work.

## Contact

Drop by my personal website [biosonix.io](https://biosonix.io) and drop a message if you would like to collaborate or need assistance with the code and development.

