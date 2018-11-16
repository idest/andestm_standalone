import numpy as np
from src.setup import input_setup
from termomecanico import termomecanico
from pyevtk.hl import imageToVTK
from src.utils import makedir

t_input, m_input = input_setup()
model, _, _ = termomecanico(t_input, m_input)

#origin=(-80, 10, -0.0498253),
#origin=(-80, 10, -0.0996506),
#spacing=(0.2, 0.2, 0.02),
makedir('vtk_files/')
imageToVTK("./vtk_files/thermomecanic_model",
        origin=(-80, 10, 0),
        spacing=(0.2, 0.2, 0.1),
        pointData={"temperature": np.asarray(model.tm.get_geotherm()),
                   "yield_strength": np.asarray(-model.mm.get_yse()[1]),
                   "geometry": np.asarray(model.gm.get_3D_geometric_model()),
                   "depth": np.asarray(model.cs.get_3D_grid()[2])})

topo = np.dstack([model.gm.get_topo().mask_irrelevant()]*1)
icd = np.dstack([model.gm.get_icd().mask_irrelevant()]*1)
moho = np.dstack([model.gm.get_moho().mask_irrelevant()]*1)
slablab = np.dstack([model.gm.get_slab_lab().mask_irrelevant()]*1)

#spacing=(0.2, 0.2, 0.02),
imageToVTK("./vtk_files/geometric_model",
        origin=(-80, 10, 0),
        spacing=(0.2, 0.2, 0.1),
        pointData={"topo": topo,
                   "icd": icd,
                   "moho": moho,
                   "slablab": slablab})
