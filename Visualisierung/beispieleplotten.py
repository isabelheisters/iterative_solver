import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Argumentwerte als 1D Arrays erzeugen
x_1d = np.linspace(0,1,301)
y_1d = np.linspace(0,1,301)

# Argumentwerte als 2D Arrays erzeugen
x_2d, y_2d = np.meshgrid(x_1d, y_1d)

# Interessante Daten erzeugen
z_2d = -(2*np.pi*4)*np.sin(x_2d*np.pi*4)*np.sin(y_2d*np.pi*4)

# Plotten
ax = plt.gca(projection='3d')
ax.plot_surface(x_2d, y_2d, z_2d, cmap=plt.get_cmap("winter"))





plt.grid()

plt.tight_layout()
plt.savefig('Visualisierung_fertig/sinableitung.png', bbox_inches='tight')