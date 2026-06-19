import matplotlib.pyplot as plt
import numpy as np

# Load data (three columns: x, y, h)
x, y, h = np.loadtxt("out.txt", unpack=True)

# Reshape to grid dimensions (33 rows, 43 columns)
X = x.reshape(33, 43)
Y = y.reshape(33, 43)
H = h.reshape(33, 43)

# Create contour plot and store the contour set
contour = plt.contour(X, Y, H)

# Add labels to the contour lines
plt.clabel(contour, inline=True, fontsize=8)

# Display the plot
plt.show()