def useScipy(ktol,ptol,bds):
	from scipy.linalg import solve
	bdNums = bds['bd1']['bdNode11']
	ktol[bdNums,bdNums] = ktol[bdNums,bdNums]*10**10
	ptol[bdNums,0] = ptol[bdNums,0]*ktol[bdNums,bdNums]
	bdNums = bds['bd1']['bdNode12']
	ktol[bdNums,bdNums] = ktol[bdNums,bdNums]*10**10
	ptol[bdNums,0] = ptol[bdNums,0]*ktol[bdNums,bdNums]
	T = solve(ktol,ptol)
	# np.savetxt("temperature.txt", np.hstack((nodes[:,0:2],T[:,0:1])))
	return T