#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import optimize_camera as optimize
import project_lidar as project

guess         = np.array([0.38, 0.026, 0.0, 977119, 210445, 145, 0.8e3])
score, params = optimize.run("lwir", num_iter=100, params=guess)
print("score  = {0}\nparams = {1}".format(score,params))

workers = 30
path    = "rawdat"
project.run_project(name="lwir", nworkers=workers, params=params, 
                    directory=path)
