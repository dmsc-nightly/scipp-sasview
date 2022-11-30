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
data.qmin = 0.0206
data.qmax = 0.294

kernel = load_model('core_shell_sphere@hardsphere')
pars = dict(radius=17,
            thickness=19,
            background=7.32,
            scale=21.84,
            sld_core=1,
            sld_shell=11,
            sld_solvent=13,
            radius_effective=37.073,
            volfraction=0.263)
model = Model(kernel, **pars)

# SET THE FITTING PARAMETERS
model.radius.range(16, 18)
model.thickness.range(18, 21)

M = Experiment(data=data, model=model)
problem = FitProblem(M)
print(f"Initial chisq {problem.chisq()}")
result = fit(problem, method='amoeba')
print(f"Final chisq {problem.chisq()}")