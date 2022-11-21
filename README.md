# scipp-sasview
Testing scipp sasview interface

### Workflow (DRAFT)
SANS2D data example in scipp/ess and pass it to sasview script for simple analysis.

1) Deploy scipp and sasview (sasmodels and sasdata to begin with) using pip and/or conda
2) Run data reduction on SANS (may require conversion of notebook to script) 
3) Store reduced data into file (requires writting NXCanSAS from scipp, there is workaround avaialable) 
4) Load file in sasview (sasdata) and do simple fitting (sasmodels)
5) Compare obtained results with expected results (check number, plots, etc)


### Known issues/things to consider
1) The worklow deploys fixed version (0.9) of ess package, otherwise micronumba install some weired dependency
2) Saving to file is currently done in ess branch. That should soon go away as https://github.com/scipp/scippneutron/issues/373#issuecomment-1306640872 is already in place
3) The idea is to checkout workflow from main, which almost always will be ahead of released version
4) SasView fit is poor and only done here for demonstration purposes
5) Need to create unit test for returned values
