import numpy as np
from src.setup import input_setup
from termomecanico import termomecanico

model = termomecanico(*input_setup())
