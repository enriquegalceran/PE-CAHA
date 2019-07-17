import numpy as np
from matplotlib.patches import Circle, Wedge, Polygon, Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

# Fixing random state for reproducibility
np.random.seed(19680801)


fig, ax = plt.subplots()

resolution = 50  # the number of vertices
N = 3
x = np.random.rand(N)
y = np.random.rand(N)
radii = 0.1*np.random.rand(N)
patches = []

# Caso 1
# patches.append(Rectangle((0, 256), 0-0+1+10, 768-256+1))
# patches.append(Rectangle((100, 256), 950-100+1, 768-256+1))
# patches.append(Rectangle((201, 513), 1900-201+1, 1536-513+1))

# Caso 2
patches.append(Rectangle((2049, 1), 2167-2049+1, 2048-1+1))
patches.append(Rectangle((1, 1), 2048-1+1, 2048-1+1))

# Caso 3
# patches.append(Rectangle((0, 1), 0-0+1+10, 2048-0+1))
# patches.append(Rectangle((501, 1), 1500-501, 2048-1+1))



ax.set_xlim(-1, 2100)
ax.set_ylim(-1, 2100)
colors = 100*np.random.rand(len(patches))
p = PatchCollection(patches, alpha=0.4)
p.set_array(np.array(colors))
ax.add_collection(p)
# fig.colorbar(p, ax=ax)



# [0,256:0,768]     	[100,256:950,768]	[201,513:1900,1536]
# [2049,1:2167,2048]	[1,1:2048,2048]  	[1,1:2048,2048]
# [0,1:0,2048]      	[501,1:1500,2048]	[501,1:1500,2048]


plt.show()
