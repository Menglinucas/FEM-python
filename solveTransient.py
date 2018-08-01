def useScipy(ktol,gtol,ptol,bds,bdParams,T0,tStart=0.,tEnd=10.e4,dt=1.0e4):
	import numpy as np
	from scipy.linalg import solve

	t = tStart
	while t <= tEnd:
		print(t)
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
		# T = np.dot(np.linalg.inv(2./3.*ktol+gtol/dt),np.dot(ptol-(1./3.*ktol-gtol/dt),T0))
		T0 = T
		t = t+dt
	# np.savetxt("temperature.txt", np.hstack((nodes[:,0:2],T[:,0:1])))
	return T