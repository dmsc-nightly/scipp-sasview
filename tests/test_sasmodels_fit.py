# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2023 dmsc-nightly contributors (https://github.com/dmsc-nightly)

from sasmodels.core import load_model
from sasmodels.bumps_model import Model, Experiment
from sasdata.dataloader.loader import Loader

from bumps.names import FitProblem
from bumps.fitters import fit

import numpy as np


def test_hardsphere_fit():

    loader = Loader()
    data = loader.load('test.nxs')

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
    result = fit(problem, method='lm')
    print(f"Final chisq {problem.chisq()}")

    #Checking if chi2 gives expected value
    np.testing.assert_almost_equal(problem.chisq(), 33.795, 3)
