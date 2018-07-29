def useScipy(ktol,ptol,bds,bdParams):
	from scipy.linalg import solve
	for key in bdParams['bd1'].keys():
		bdNums = bds['bd1'][key+'Node']
		ktol[bdNums,bdNums] = ktol[bdNums,bdNums]*10**10
		ptol[bdNums,0] = ptol[bdNums,0]*ktol[bdNums,bdNums]
	T = solve(ktol,ptol)
	# np.savetxt("temperature.txt", np.hstack((nodes[:,0:2],T[:,0:1])))
	return T