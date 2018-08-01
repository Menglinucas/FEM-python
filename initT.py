def initToBeZero(nodes,bds,bdParams):
	import numpy as np
	T0 = np.zeros((len(nodes),1))
	# add the 1st boundary
	if 'bd1' in bdParams.keys():
		for key in bdParams['bd1'].keys():
			if key == 'bdpts':
				# note the dimension of T0 and bds[..]
				T0[bds['bd1'][key+'Node']] = np.reshape(bds['bd1'][key+'Value'],(len(T0[bds['bd1'][key+'Node']]),1))
			else:
				T0[bds['bd1'][key+'Node']] = bdParams['bd1'][key]
	return T0