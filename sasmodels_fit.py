from sasmodels.core import load_model
from sasmodels.bumps_model import Model, Experiment
#from sasmodels.data import load_data
from sasdata.dataloader.loader import Loader

from bumps.names import *
from bumps.fitters import fit
from bumps.formatnum import format_uncertainty


loader = Loader()
data = loader.load('test.h5')
#data = load_data('test_data.txt')

data = data[0]
data.qmin = min(data.x)
data.qmax = max(data.x)

kernel = load_model('cylinder')
pars = dict(radius=35,
            length=350,
            background=0.0,
            scale=1.0,
            sld=4.0,
            sld_solvent=1.0)
model = Model(kernel, **pars)

# SET THE FITTING PARAMETERS
model.radius.range(1, 50)
model.length.range(1, 500)

M = Experiment(data=data, model=model)
problem = FitProblem(M)
print(f"Initial chisq {problem.chisq()}")
result = fit(problem, method='amoeba')
print(f"Final chisq {problem.chisq()}")