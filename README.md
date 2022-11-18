# scipp-sasview
Testing scipp sasview interface

### Workflow (DRAFT)
SANS2D data example in scipp/ess and pass it to sasview script for simple analysis.

1) Deploy scipp and sasview (sasmodels and sasdata to begin with) using pip
2) Run data reduction on SANS (may require conversion of notebook to script) 
3) Store reduced data into file (requires writting NXCanSAS from scipp, there is workaround avaialable) 
4) Load file in sasview (sasdata) and do simple fitting (sasmodels)
5) Compare obtained results with expected results (check number, plots, etc)
