import sys
import numpy as np
from matplotlib import pyplot

data = np.loadtxt(sys.argv[1]);

color_map = ["red", "green", "blue"]
labels = data[:, 2]
label_colors = [ color_map[int(label)] for label in labels ]

pyplot.scatter(data[:,0], data[:,1], c=label_colors)
pyplot.show()
