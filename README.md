![scipp-sasview](https://github.com/dmsc-nightly/scipp-sasview/actions/workflows/nightly.yml/badge.svg?branch=main)

# scipp-sasview
Testing the Scipp-Sasview interface

### Workflow

Use the SANS2D [I(Q) workflow](https://scipp.github.io/ess/instruments/loki/sans2d_to_I_of_Q.html) from scipp/ess and feed the results to a Sasview script for simple analysis.

1. Install Scipp and Sasview (sasmodels and sasdata to begin with) using `pip` and `conda`
1. Convert the Jupyter notebook that contains the reduction workflow to an executable script.
1. Run the SANS data reduction and store reduced data into file (requires writing NXCanSAS from Scipp).
1. Load file in Sasview (sasdata) and do simple fitting (sasmodels).
1. Compare obtained results with expected results (check chisquared, plots, etc).
