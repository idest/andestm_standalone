import numpy as np
import matplotlib.pyplot as plt
from src.colormaps import categorical_cmap

c1 = categorical_cmap(4, 3, cmap="tab10")
plt.scatter(np.arange(4*3),np.ones(4*3)+1, c=np.arange(4*3), s=180, cmap=c1)

c2 = categorical_cmap(2, 5, cmap="tab10")
plt.scatter(np.arange(10),np.ones(10), c=np.arange(10), s=180, cmap=c2)

c3 = categorical_cmap(5, 4, cmap="tab10")
plt.scatter(np.arange(20),np.ones(20)-1, c=np.arange(20), s=180, cmap=c3)    

c4 = categorical_cmap(3, [3,5,3])
plt.scatter(np.arange(11), np.ones(11)-2, c=np.arange(11), s=180, cmap=c4)

plt.margins(y=0.3)
plt.xticks([])
plt.yticks([0,1,2],["(5, 4)", "(2, 5)", "(4, 3)"])
plt.show()
