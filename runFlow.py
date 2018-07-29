def main():
	import preProcess as prep
	import buildStiffMatrix as bsm
	import solveEquation as solEq
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
	
	# materials
	matParams = {'mat1':[3.,3000.,10000.,1e-4,1e-4,1.2e-6],
				'mat2':[10.,3000.,10000.,0,0,0]}
	# boundaries
	bdParams = {'bd1':{'bd11':10.,
						'bd12':100.},
				'bd2':{'bd21':0.}}
	# bdParams = {'bd1':{'bd11':10.,
						# 'bd12':100.},
				# 'bd2':{'bd21':0.},
				# 'bd3':{'bd31':[3.,0.]}}

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
	ktol,ptol = bsm.tolStiff(nodes,mats,bds,bdParams)

	##########################################################################
	######################## 3. solve linear equations #######################
	##########################################################################
	T = solEq.useScipy(ktol,ptol,bds,bdParams)

	##########################################################################
	############################# 4. postprocess #############################
	########################################################################## 
	postp.drawHeat(x=nodes[:,0],y=nodes[:,1],z=T)
	
if __name__ == '__main__':
	main()