def drawHeatByLinear(nodes,values):
	from scipy.interpolate import Rbf
	from scipy.interpolate import griddata
	from matplotlib import cm
	import matplotlib.pyplot as plt
	import numpy as np
	fig = plt.figure()
	x = nodes[:,0]; y = nodes[:,1]; z = values
	z = np.reshape(z,(len(z)))
	gridx, gridy = np.mgrid[np.min(x):np.max(x):200j, np.min(y):np.max(y):200j]
	gridz = griddata(nodes[:,0:2],z,(gridx,gridy),method='linear')
	plt.imshow(gridz.T,extent=(np.min(x),np.max(x),np.min(y),np.max(y)),origin='lower',cmap=cm.jet)
	plt.title('FEM test')
	plt.colorbar()
	plt.show()
	
def drawHeatByRbf(nodes,values):
	from scipy.interpolate import Rbf
	from matplotlib import cm
	import matplotlib.pyplot as plt
	import numpy as np
	fig = plt.figure()
	x = nodes[:,0]; y = nodes[:,1]; z = values
	gridx = np.linspace(np.min(x),np.max(x),200)
	gridy = np.linspace(np.min(y),np.max(y),200)
	XI,YI = np.meshgrid(gridx,gridy)
	rbf = Rbf(x,y,z,epsilon=2)
	ZI = rbf(XI,YI)
	plt.pcolor(XI,YI,ZI,cmap=cm.jet)
	plt.title('FEM test')
	plt.xlim(np.min(x),np.max(x))
	plt.ylim(np.min(y),np.max(y))
	plt.colorbar()
	plt.show()