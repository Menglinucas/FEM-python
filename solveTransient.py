def useScipy(nodes,ktol,gtol,ptol,bds,bdParams,T0,tStart,tEnd,dt):
	import numpy as np
	from scipy.linalg import solve
	from scipy.interpolate import Rbf
	from scipy.interpolate import griddata
	from matplotlib import cm
	import matplotlib.pyplot as plt
	import numpy as np

	plt.close()
	fig = plt.figure()
	x = nodes[:,0]; y = nodes[:,1]
	gridx, gridy = np.mgrid[np.min(x):np.max(x):200j, np.min(y):np.max(y):200j]
	plt.ion()
	
	t = tStart
	while t <= tEnd:
		print(t/3.1536e13)
		leftArray = 2./3.*ktol+gtol/dt
		rightArray = ptol-np.dot((1./3.*ktol-gtol/dt),T0)
		# add the 1st boundary
		if 'bd1' in bdParams.keys():
			for key in bdParams['bd1'].keys():
				if key == 'bdpts':
					bdNums = bds['bd1'][key+'Node']
					leftArray[bdNums,bdNums] = leftArray[bdNums,bdNums]*10**10
					rightArray[bdNums,0] = bds['bd1'][key+'Value']*leftArray[bdNums,bdNums]
				else:
					bdNums = bds['bd1'][key+'Node']
					leftArray[bdNums,bdNums] = leftArray[bdNums,bdNums]*10**10
					rightArray[bdNums,0] = bdParams['bd1'][key]*leftArray[bdNums,bdNums]
		T = solve(leftArray,rightArray)
		# draw the heat map
		z = T
		z = np.reshape(z,(len(z)))
		gridz = griddata(nodes[:,0:2],z,(gridx,gridy),method='linear')
		plt.title(str(t/3.1536e13)+'Ma')
		plt.imshow(gridz.T,extent=(np.min(x),np.max(x),np.min(y),np.max(y)),origin='lower',cmap=cm.jet)
		plt.pause(0.01)
		plt.cla()
		# T = np.dot(np.linalg.inv(2./3.*ktol+gtol/dt),np.dot(ptol-(1./3.*ktol-gtol/dt),T0))
		T0 = T
		t = t+dt
	# np.savetxt("temperature.txt", np.hstack((nodes[:,0:2],T[:,0:1])))
	return T