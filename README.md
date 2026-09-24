# Client Exchange Dynamics

MATLAB code for the manuscript  
**“The Exchange Dynamics of Client Molecules in Biomolecular Condensates”**  
by Ross Kliegman, Vladimir Grigorev, and Yaojun Zhang.

## Files

- `Data_Main.m` — numerically solves the two-state reaction-diffusion model and saves the results in `PDEresults.mat`.
- `Figure2_plot.m` — generates Fig. 2 using `PDEresults.mat`.
- `Figure3_plot.m` — generates Fig. 3 using `PDEresults.mat`.

## Usage

First run

```matlab
Data_Main
```

to generate `PDEresults.mat`.

> **Note:** The generated `PDEresults.mat` file is approximately 5 GB in size and may take substantial time and memory to generate.

Then run

```matlab
Figure2_plot
```

or

```matlab
Figure3_plot
```

to reproduce Figs. 2 and 3, respectively.

## Requirements

MATLAB with the Curve Fitting Toolbox.

## License

This code is distributed under the MIT License.
