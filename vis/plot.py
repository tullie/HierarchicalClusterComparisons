import sys
import random
import numpy as np
from matplotlib import pyplot

data = np.loadtxt(sys.argv[1]);

color_map = ["red", "green", "blue", "yellow", "black", "magenta", "cyan", "#fb09da", "#c3903b", "#bd0cfb", "#c22014", "#e3abc5", '#e3abc5', '#c22014', '#2deead', '#717157']
labels = data[:, 2]

#label_colors = []
#for i in range(0, len(labels)):
#    if i < len(color_map):
#        label_colors.append(color_map[int(labels[i])])
#    else:
#        r = lambda: random.randint(0,255)
#        hex = '#%02X%02X%02X' % (r(),r(),r())
#        label_colors.append(hex)

pyplot.scatter(data[:,0], data[:,1], c=labels)
pyplot.show()
