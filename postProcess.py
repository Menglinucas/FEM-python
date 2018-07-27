def drawHeat(x,y,z):
	from scipy.interpolate import Rbf
	from matplotlib import cm
	import matplotlib as plt
	import matplotlib.pyplot as plt
	import numpy as np
	gridx = np.linspace(0,10,200)
	gridy = np.linspace(0,10,200)
	XI,YI = np.meshgrid(gridx,gridy)
	rbf = Rbf(x,y,z,epsilon=2)
	ZI = rbf(XI,YI)
	plt.subplot(1,1,1)
	plt.pcolor(XI,YI,ZI,cmap=cm.jet)
	plt.title('FEM test')
	plt.xlim(0,10)
	plt.ylim(0.10)
	plt.colorbar()
	plt.show()