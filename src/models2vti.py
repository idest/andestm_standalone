import numpy as np
from setup import input_setup
from termomecanico import termomecanico
from pyevtk.hl import imageToVTK
from utils import makedir

t_input, m_input = input_setup()
model, _, _ = termomecanico(t_input, m_input)

#origin=(-80, 10, -0.0498253),
makedir('vtk_files/')
imageToVTK("./vtk_files/thermomecanic_model",
        origin=(-80, 10, -0.0996506),
        spacing=(0.2, 0.2, 0.02), 
        pointData={"temperature": model.tm.get_geotherm().array,
                   "yield_strength": model.mm.get_yse()[0].array})

topo = np.dstack([model.gm.get_topo().mask_irrelevant(nan_fill=True)]*1)
icd = np.dstack([model.gm.get_icd().mask_irrelevant(nan_fill=True)]*1)
moho = np.dstack([model.gm.get_moho().mask_irrelevant(nan_fill=True)]*1)
slablab = np.dstack([model.gm.get_slab_lab().mask_irrelevant(nan_fill=True)]*1)

imageToVTK("./vtk_files/geometric_model",
        origin=(-80, 10, 0),
        spacing=(0.2, 0.2, 0.02),
        pointData={"topo": topo,
                   "icd": icd,
                   "moho": moho,
                   "slablab": slablab})
