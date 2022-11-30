from sasmodels.core import load_model
from sasmodels.bumps_model import Model, Experiment
from sasdata.dataloader.loader import Loader

from bumps.names import *
from bumps.fitters import fit
#from bumps.formatnum import format_uncertainty


loader = Loader()
data = loader.load('test.nxs')
#data = load_data('test_data.txt')

data = data[0]
data.qmin = 0.0229
data.qmax = 0.151

kernel = load_model('sphere@hardsphere')
pars = dict(radius=29,
            background=18.0,
            scale=2.96,
            sld=3.83,
            sld_solvent=14.73,
            radius_effective=36.72,
            volfraction=0.258)
model = Model(kernel, **pars)

# SET THE FITTING PARAMETERS
model.radius.range(25, 35)

M = Experiment(data=data, model=model)
problem = FitProblem(M)
print(f"Initial chisq {problem.chisq()}")
result = fit(problem, method='amoeba')
print(f"Final chisq {problem.chisq()}")