# Client Exchange Dynamics

MATLAB code for the manuscript  
**“The Exchange Dynamics of Client Molecules in Biomolecular Condensates”**  
by Ross Kliegman, Vladimir Grigorev, and Yaojun Zhang.

## Files

- `Data_Main.m` — numerically solves the full two-state reaction-diffusion model for bound and unbound clients and saves the results in `PDEresults.mat`.
- `Data_Main_SingleState.m` — numerically solves the effective single-state model and saves the results in `PDEresults_SingleState.mat`.
- `Figure2_plot.m` — generates Fig. 2 using `PDEresults.mat`.
- `Figure3_plot.m` — generates Fig. 3 using `PDEresults.mat`.
- `Figure4_plot.m` — generates Fig. 4 by comparing the full two-state model in `PDEresults.mat` with the effective single-state model in `PDEresults_SingleState.mat`.

## Usage

First run

```matlab
Data_Main
```

to generate `PDEresults.mat`.

To generate the effective single-state results used in Fig. 4, also run

```matlab
Data_Main_SingleState
```

to generate `PDEresults_SingleState.mat`.

> **Note:** The generated MATLAB data files are large. In particular, `PDEresults.mat` is approximately 5 GB and may take substantial time and memory to generate.

Then run

```matlab
Figure2_plot
```

or

```matlab
Figure3_plot
```

to reproduce Figs. 2 and 3, respectively.

After both `PDEresults.mat` and `PDEresults_SingleState.mat` have been generated, run

```matlab
Figure4_plot
```

to reproduce Fig. 4.

## Requirements

MATLAB with the Curve Fitting Toolbox.

## License

This code is distributed under the MIT License.
