import numpy as np
# tri: [node1,node2,node3]; nodes: [x,y,z]
def triShapeFunc(nodes,tri):
	x1 = nodes[tri[0]][0]; y1 = nodes[tri[0]][1]
	x2 = nodes[tri[1]][0]; y2 = nodes[tri[1]][1]
	x3 = nodes[tri[2]][0]; y3 = nodes[tri[2]][1]
	a1 = y2-y3; b1 = x3-x2; c1 = x2*y3-x3*y2
	a2 = y3-y1; b2 = x1-x3; c2 = x3*y1-x1*y3
	a3 = y1-y2; b3 = x2-x1; c3 = x1*y2-x2*y1
	delta = (a1*b2-a2*b1)/2.0
	# n1 = [a1/area/2.0, b1/area/2.0, c1/area/2.0] # n1 = L1 = (a1*x+b1*y+c1)/2/area
	# n2 = [a2/area/2.0, b2/area/2.0, c2/area/2.0] # n2 = L2 = (a2*x+b2*y+c2)/2/area
	# n3 = [a3/area/2.0, b3/area/2.0, c3/area/2.0] # n3 = L3 = (a3*x+b3*y+c3)/2/area
	# in postprocess: u = n1u1 + n2u2 + n3u3
	return x1,x2,x3,y1,y2,y3,a1,a2,a3,b1,b2,b3,c1,c2,c3,delta
	
# params1: kappa,miu,miuW,vx,vy,Q; alpha, beita
def boundParams(tri, bds, bdParams):
	# whether containing a second or third boundary
	isB = []; notB = []; eleAlpha = 0.; eleBeita = 0.
	anyTwo = [[0,1],[0,2],[1,2]]
	for keyDict in bdParams.keys():
		for key in bdParams[keyDict]:
			if key == 'bdpts': continue
			if set(tri[anyTwo[0]]).issubset(bds[keyDict][key+'Node']):
				isB = anyTwo[0]; eleAlpha, eleBeita = bds[keyDict][key+'Params']
			elif set(tri[anyTwo[1]]).issubset(bds[keyDict][key+'Node']):
				isB = anyTwo[1]; eleAlpha, eleBeita = bds[keyDict][key+'Params']
			elif set(tri[anyTwo[2]]).issubset(bds[keyDict][key+'Node']):
				isB = anyTwo[2]; eleAlpha, eleBeita = bds[keyDict][key+'Params']
	notB = list(set([0,1,2])-set(isB))
	return isB, notB, eleAlpha, eleBeita
	
# elementary stiffness
def eleStiff(nodes,tri,bds,bdParams,eleK,eleMiu,eleMiuW,eleVx,eleVy,Q): 
	x1,x2,x3,y1,y2,y3,a1,a2,a3,b1,b2,b3,c1,c2,c3,delta = triShapeFunc(nodes,tri)
	x = [x1,x2,x3]; y = [y1,y2,y3]
	isB, notB, eleAlpha, eleBeita = boundParams(tri, bds, bdParams)
	# K
	ke1 = eleK/4./delta*np.array([[a1*a1+b1*b1, a1*a2+b1*b2, a1*a3+b1*b3],
									[a2*a1+b2*b1, a2*a2+b2*b2, a2*a3+b2*b3],
									[a3*a1+b3*b1, a3*a2+b3*b2, a3*a3+b3*b3]])
	ke2 = eleMiuW/6.*np.array([[eleVx*a1+eleVy*b1, eleVx*a1+eleVy*b2, eleVx*a1+eleVy*b3],
							[eleVx*a2+eleVy*b1, eleVx*a2+eleVy*b2, eleVx*a2+eleVy*b3],
							[eleVx*a3+eleVy*b1, eleVx*a3+eleVy*b2, eleVx*a3+eleVy*b3]])	
	ke3 = np.zeros((3,3))	
	ge = eleMiu*delta/12.*np.array([[2,1,1],
								[1,2,1],
								[1,1,2]])
	pe1 = Q*delta/3.*np.ones((3,1))
	pe2 = np.zeros((3,1))
	# deal with the boundary
	if len(isB) == 2:
		ke3 = eleAlpha/6.*np.sqrt((x[isB[1]]-x[isB[0]])**2+(y[isB[1]]-y[isB[0]])**2)*np.array([[2,1,1],
																								[1,2,1],
																								[1,1,2]])
		ke3[notB[0],:] = 0; ke3[:,notB[0]] = 0
		pe2 = eleBeita/2.*np.sqrt((x[isB[1]]-x[isB[0]])**2+(y[isB[1]]-y[isB[0]])**2)*np.array([[1],
																								[1],
																								[1]])
		pe2[notB[0],0] = 0
	# preliminary sum
	ke = ke1+ke2+ke3
	pe = pe1+pe2
	return ke,ge,pe

# total stiffness
def tolStiff(nodes,mats,bds,bdParams):
	length = len(nodes)
	ktol = np.zeros((length,length))
	gtol = np.zeros((length,length))
	ptol = np.zeros((length,1))
	for mat in mats.keys():
		for tri in mats[mat]['tris']:
			kappa,miu,miuW,vx,vy,Q = mats[mat]['params']
			ke,ge,pe = eleStiff(nodes,tri,bds,bdParams,kappa,miu,miuW,vx,vy,Q)
			ktol[tri[0],tri[0]] = ktol[tri[0],tri[0]]+ke[0,0]
			ktol[tri[0],tri[1]] = ktol[tri[0],tri[1]]+ke[0,1]
			ktol[tri[0],tri[2]] = ktol[tri[0],tri[2]]+ke[0,2]
			ktol[tri[1],tri[0]] = ktol[tri[1],tri[0]]+ke[1,0]
			ktol[tri[1],tri[1]] = ktol[tri[1],tri[1]]+ke[1,1]
			ktol[tri[1],tri[2]] = ktol[tri[1],tri[2]]+ke[1,2]
			ktol[tri[2],tri[0]] = ktol[tri[2],tri[0]]+ke[2,0]
			ktol[tri[2],tri[1]] = ktol[tri[2],tri[1]]+ke[2,1]
			ktol[tri[2],tri[2]] = ktol[tri[2],tri[2]]+ke[2,2]
			gtol[tri[0],tri[0]] = gtol[tri[0],tri[0]]+ke[0,0]
			gtol[tri[0],tri[1]] = gtol[tri[0],tri[1]]+ke[0,1]
			gtol[tri[0],tri[2]] = gtol[tri[0],tri[2]]+ke[0,2]
			gtol[tri[1],tri[0]] = gtol[tri[1],tri[0]]+ke[1,0]
			gtol[tri[1],tri[1]] = gtol[tri[1],tri[1]]+ke[1,1]
			gtol[tri[1],tri[2]] = gtol[tri[1],tri[2]]+ke[1,2]
			gtol[tri[2],tri[0]] = gtol[tri[2],tri[0]]+ke[2,0]
			gtol[tri[2],tri[1]] = gtol[tri[2],tri[1]]+ke[2,1]
			gtol[tri[2],tri[2]] = gtol[tri[2],tri[2]]+ke[2,2]
			ptol[tri[0],0] = ptol[tri[0],0]+pe[0,0]
			ptol[tri[1],0] = ptol[tri[1],0]+pe[1,0]
			ptol[tri[2],0] = ptol[tri[2],0]+pe[2,0]
	return ktol,gtol,ptol