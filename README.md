# Mechanistic Modeling of CRISPRi/CRISPRa

This repository contains supplementary information for the paper: 
> Scott, Helen, Alessandro Occhialini, Scott C. Lenaghan, and Jacob Beal. 2024.
"Simulations Predict Stronger CRISPRi Transcriptional Repression in Plants for
Identical than Heterogeneous Target Sites." bioRxiv.
https://doi.org/10.1101/2024.04.22.590637

If you make use of any of the contents of this repository, please cite this paper.

## Abstract
Plant synthetic biologists have been working to adapt the CRISPRa and CRISPRi
promoter regulation methods for applications such as improving crops or
installing other valuable pathways. With other organisms, strong
transcriptional control has typically required multiple gRNA target sites,
which poses a critical engineering choice between heterogeneous sites, which
allow each gRNA to target existing locations in a promoter, and identical
sites, which typically require modification of the promoter. Here we
investigate the consequences of this choice for CRISPRi plant promoter
regulation via simulation-based analysis, using model parameters based on
single gRNA regulation and constitutive promoters in _N. benthamiana_.
Using models of 2 to 6 gRNA target sites to compare heterogeneous versus
identical sites for tunability, sensitivity to parameter values, and
sensitivity to cell-to-cell variation, we find that identical gRNA target
sites are predicted to yield far more effective transcriptional repression than
heterogeneous sites.

## Contents of Repository
- Python code that generates SBOL models each potential promoter architecture.
- Python code to export SBOL models into ODE equations in LaTeX and executable ODE simulation models in MATLAB.
- Simulation runners for experimentation with the models.
- Results and figures generated from these simulations.

## Running the code

### How to Make a Circuit Model

The routines for building SBOL circuit components are in `sbol/builders.py`.
Import the module and run via a script like `sbol/make_sbol_models.py`.
The resulting circuit models we generated are saved in `sbol/gRNA_models.nt`.

### Generating LaTeX Equations

The routines for generating LaTeX from SBOL circuits are in `sbol/latex_generation.py`.
Import the module and run via a script like `sbol/sbol_to_latex.py`.
The generated LaTeX equations can be embedded in a LaTeX document like `equations/view_generated_equations.tex`.
The resulting LaTeX equations we generated are saved in `equations/generated_equations.tex` and are viewable in `equations/view_generated_equations.pdf`.

### Generating MATLAB Code

The routines for generating MATLAB code from SBOL circuits are in `sbol/matlab_generation.py`.
Import the module and run via a script like `sbol/sbol_to_matlab.py`.
The resulting MATLAB models we generated are saved in the `models/` directory.
