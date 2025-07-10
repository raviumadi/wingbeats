# Wingbeat-Call Asynchrony Simulation and Theory

This repository provides the full simulation code, analysis scripts, and supporting functions used in the study of **Temporal Precision** built on the **Responsivity Framework** in echolocating bats. It contains modules to explore the temporal dynamics of sonar call timing, wingbeat coordination, and sensory feedback loops in echolocating animals.

This code supports the following manuscripts:

1. [*Biosonar Responsivity Sets the Stage for the Terminal Buzz* ](https://www.biorxiv.org/content/10.1101/2025.06.16.659925v1) 
2. [*Temporal Precision Necessitates Wingbeat-Call Asynchrony in Actively Echolocating Bats*](https://www.biorxiv.org/content/10.1101/2025.06.18.660328v1)  
3. [*Swarm Cohesion in Bats Emerges from Stable Temporal Loops*](https://www.biorxiv.org/content/10.1101/2025.07.05.663265v1)

---

## 📁 Repository Structure

### `/src`  
Core functions for signal generation, propagation modelling, and acoustic analysis:

- `simulateEcholocationWings.m` – primary simulation engine for sonar–wingbeat coupling  
- `generateVirtualBatCall.m` – generates frequency-modulated bat calls  
- `airAttenuationFilter.m`, `calculateTargetStrength.m`, `calculateGeometricAttenuation.m`, etc. – physical modelling of sound propagation  
- `formatLatex*.m` – visual formatting for plots with LaTeX labels  

### `/mat`  
Simulation entry points and experiments scripts:

- `wingbeat_sim.m`, `wingbeat_sim_loop.m` – test cases and multi-run simulations for wingbeat–call synchrony  
- `sim_rand_var_run.m` – Monte Carlo simulations with randomly sampled biological parameters  
- `ssg_demo.m` – sonar sound group simulation demo  
- `supplementary_analysis.m` – statistical and graphical analyses for paper figures  
- `temporal_integrity.m` – analysis of call–echo timing fidelity  

### `/results`

Folder to store figures, `.mat` files and output summaries from simulations and data processing.  

### `/data`  
Contains field recordings and/or preprocessed metrics, where applicable.

### `/doc`  
Documentation, paper manuscripts, and supporting LaTeX or exportable figures.

---

## 🧪 Getting Started

Ensure you have MATLAB R2020b or later. Run the following script to add all folders to the path:

```matlab
init
```



## **📜 License**

This project is licensed under the [Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)](https://creativecommons.org/licenses/by-nc/4.0/).

You are free to share and adapt the material for non-commercial use with appropriate credit.

## **🧠 Citation**

If you use this codebase in part or in full, please cite the corresponding manuscripts or contact the author for further information.

## Other Publications

Check my [ResearchGate](https://www.researchgate.net/profile/Ravi-Umadi-3) or [ORCID](https://orcid.org/0000-0003-3867-1769) for a full list of my research work. 

## Contact

Drop by my personal website [biosonix.io](https://biosonix.io) and drop a message if you would like to collaborate or need assistance with the code and development. 

