import sys
import numpy as np
from matplotlib import pyplot

data = np.loadtxt(sys.argv[1]);
pyplot.scatter(data[:,0], data[:,1], c=data[:, 2])
pyplot.show()
