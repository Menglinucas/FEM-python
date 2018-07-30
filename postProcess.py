def drawHeatByLinear(nodes,values):
	from scipy.interpolate import Rbf
	from scipy.interpolate import griddata
	from matplotlib import cm
	import matplotlib.pyplot as plt
	import numpy as np
	x = nodes[:,0]; y = nodes[:,1]; z = values
	z = np.reshape(z,(len(z)))
	gridx, gridy = np.mgrid[np.min(x):np.max(x):200j, np.min(y):np.max(y):200j]
	gridz = griddata(nodes[:,0:2],z,(gridx,gridy),method='linear')
	plt.subplot(1,1,1)
	plt.imshow(gridz.T,extent=(np.min(x),np.max(x),np.min(y),np.max(y)),origin='lower',cmap=cm.jet)
	plt.colorbar()
	plt.show()
	
def drawHeatByRbf(nodes,values):
	from scipy.interpolate import Rbf
	from matplotlib import cm
	import matplotlib as plt
	import matplotlib.pyplot as plt
	import numpy as np
	x = nodes[:,0]; y = nodes[:,1]; z = values
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