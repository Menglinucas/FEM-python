def main(trans=False,tStart=0.,tEnd=1.0e14,dt=1.0e13):
	import numpy as np
	import preProcess as prep
	import buildStiffMatrix as bsm
	import solveStatic as solS
	import initT
	import solveTransient as solT
	import postProcess as postp

	##########################################################################
	############################## 1. preprocess #############################
	##########################################################################
	# (1) set parameters
	# kappa --- conductivity
	# miu --- specific heat capacity of the media
	# miuW --- specific heat capacity of the fluid
	# alpha --- heat exchange coefficient
	# Ts --- the temperature of exchange heat source
	# Q --- heat source
	# vx --- direct-x velocity
	# vy --- direct-y velocity
	# qhf --- heat flow
	# beita --- a parameter
	
	# materials, [kappa, miu, miuW, vx, vy, Q]
	matParams = {'mat1':[2.7,3.e6,0.,1./100./24./3600.,1./100./24./3600.,1.e-4],
				'mat2':[3.2,3.e6,5.e6,0.,0.,0.]}
	# boundaries
	bdParams = {'bd1':{'bd11':1200.,
						'bd12':0.,    # T
						'bdpts':[[2.5e3,2.5e3,100],[7.5e3,7.5e3,900]]},   # [[x,y,T],...]
				'bd2':{'bd21':0.1},         # q
				'bd3':{'bd31':[3.,1200.]}}  # [alpha, beita]

	# (2) generate mesh (command or interface)
	mesh = prep.setModelByCommand(meshPath='theMesh/theMesh.msh')
	# mesh = pp.setModelByInterface()
	
	# (3) get nodes, triangles, boundaries
	nodes = prep.getNodes(mesh)     
	bds = prep.getBoundaries(mesh,bdParams) 
	mats = prep.getSurfsWithTrisAndParams(mesh,matParams)

	##########################################################################
	######################## 2. build stiffness matrix #######################
	##########################################################################
	ktol,gtol,ptol = bsm.tolStiff(nodes,mats,bds,bdParams)
	
	##########################################################################
	######################## 3. solve linear equations #######################
	##########################################################################
	if trans == False:
		# static
		T = solS.useScipy(ktol,ptol,bds,bdParams)
	else:
		# transient
		T0 = initT.initToBeZero(nodes,bds,bdParams)
		T = solT.useScipy(ktol,gtol,ptol,bds,bdParams,T0,tStart=tStart,tEnd=tEnd,dt=dt)
	
	##########################################################################
	############################# 4. postprocess #############################
	########################################################################## 
	postp.drawHeatByRbf(nodes=nodes,values=T)
	
if __name__ == '__main__':
	main()