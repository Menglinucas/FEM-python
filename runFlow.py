def main():
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
	matParams = {'mat1':[2.,3000.,0.,0.,0.,1.e-3],
				'mat2':[3.,3000.,10000.,0.e-4,0.e-4,0.]}
	# boundaries
	bdParams = {'bd1':{'bd11':0.,
						'bd12':100.,    # T
						'bdpts':[[2.5,2.5,90],[7.5,7.5,10]]},   # [[x,y,T],...]
				'bd2':{'bd21':10.},         # q
				'bd3':{'bd31':[0.,0.]}}  # [alpha, beita]

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
	# static
	# T = solS.useScipy(ktol,ptol,bds,bdParams)
	
	# transient
	T0 = initT.initToBeZero(nodes)
	T = solT.useScipy(ktol,gtol,ptol,bds,bdParams,T0,tStart=0.,tEnd=1.e2,dt=1.0e1)
	
	##########################################################################
	############################# 4. postprocess #############################
	########################################################################## 
	postp.drawHeatByRbf(nodes=nodes,values=T)
	
if __name__ == '__main__':
	main()