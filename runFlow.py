def main():
	import preProcess as prep
	import buildStiffMatrix as bsm
	import solveEquation as solEq
	import postProcess as postp

	##########################################################################
	############################## 1. preprocess #############################
	##########################################################################
	# (1) set parameters
	kappa=3.0    # conductivity
	miu=3000.    # specific heat capacity of the media
	miuW=10000.  # specific heat capacity of the fluid
	alpha=3.     # heat exchange coefficient
	Ts=0.        # the temperature of exchange heat source
	Q=0.         # heat source
	vx=0.        # direct-x velocity
	vy=0.        # direct-y velocity
	qhf=0.6      # heat flow
	beita=0.     # a parameter

	# (2) generate mesh (command or interface)
	mesh = prep.setModelByCommand(meshPath='theMesh/theMesh.msh')
	# mesh = pp.setModelByInterface()
	
	# (3) get nodes, triangles, boundaries
	nodes = prep.getNodes(mesh)    # global variable
	tris = prep.getTriangles(mesh) # global variable
	bds = prep.getBoundaries(mesh,alpha,Ts) # global variable

	##########################################################################
	######################## 2. build stiffness matrix #######################
	##########################################################################
	ktol,ptol = bsm.tolStiff(nodes,tris,bds,kappa,alpha,beita,miu,miuW,vx,vy,Q)

	##########################################################################
	######################## 3. solve linear equations #######################
	##########################################################################
	T = solEq.useScipy(ktol,ptol,bds)

	##########################################################################
	############################# 4. postprocess #############################
	########################################################################## 
	postp.drawHeat(x=nodes[:,0],y=nodes[:,1],z=T)
	
if __name__ == '__main__':
	main()