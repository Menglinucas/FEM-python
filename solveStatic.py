def useScipy(ktol,ptol,bds,bdParams):
	from scipy.linalg import solve
	import numpy as np
	# add the 1st boundary
	if 'bd1' in bdParams.keys():
		for key in bdParams['bd1'].keys():
			if key == 'bdpts':
				bdNums = bds['bd1'][key+'Node']
				ktol[bdNums,bdNums] = ktol[bdNums,bdNums]*10**10
				ptol[bdNums,0] = bds['bd1'][key+'Value']*ktol[bdNums,bdNums]
			else:
				bdNums = bds['bd1'][key+'Node']
				ktol[bdNums,bdNums] = ktol[bdNums,bdNums]*10**10
				ptol[bdNums,0] = bdParams['bd1'][key]*ktol[bdNums,bdNums]
	T = solve(ktol,ptol)
	# T = np.dot(np.linalg.inv(ktol),ptol)
	# np.savetxt("temperature.txt", np.hstack((nodes[:,0:2],T[:,0:1])))
	return T
	
	
	