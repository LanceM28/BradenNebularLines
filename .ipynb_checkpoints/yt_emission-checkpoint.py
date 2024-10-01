# importing packages
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

import yt
# Load the dataset.
ds = yt.load("/Users/bnowicki/Documents/Research/Ricotti/output_00273")

print(dir(ds.fields.gas))

cell_fields = [
    "Density",
    "x-velocity",
    "y-velocity",
    "z-velocity",
    "Pressure",
    "Metallicity",
    # "dark_matter_density",
    "xHI",
    "xHII",
    "xHeII",
    "xHeIII",
]
epf = [
    ("particle_family", "b"),
    ("particle_tag", "b"),
    ("particle_birth_epoch", "d"),
    ("particle_metallicity", "d"),
]



f1 = "/Users/bnowicki/Documents/Research/Ricotti/output_00273"

ds = yt.load(f1, fields=cell_fields, extra_particle_fields=epf)

p = yt.ProjectionPlot(ds, "z", ("gas", "number_density"), width=0.0001,
                      weight_field=("gas", "number_density"),
                      buff_size=(1000, 1000),
                      center=[0.49118094, 0.49275361, 0.49473726])

p.show()
