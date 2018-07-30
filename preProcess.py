# gmsh path
# gmshPath = 'D:/software/gmsh-3.0.6-Windows64/gmsh.exe'

# set the model
def setModelByInterface():
	import os
	# os.system(gmshPath)
	os.system('gmsh')
	
def loadMesh(meshPath='theMesh/theMesh.msh'):
	import meshio
	mesh = meshio.read(meshPath)
	return mesh

def setModelByCommand(meshPath='theMesh/theMesh.msh'):
	import pygmsh as pg
	import meshio
	import numpy as np
	geom = pg.built_in.Geometry()
	# add point
	pts = {}
	pts[1] = geom.add_point([0,0,0],lcar=1)
	pts[2] = geom.add_point([10,0,0],lcar=1)
	pts[3] = geom.add_point([10,10,0],lcar=1)
	pts[4] = geom.add_point([0,10,0],lcar=1)
	pts[5] = geom.add_point([10,5,0],lcar=1)
	pts[6] = geom.add_point([0,5,0],lcar=1)
	# add line
	lines = {}
	lines[1] = geom.add_line(pts[1],pts[2])
	lines[2] = geom.add_line(pts[2],pts[5])
	lines[3] = geom.add_line(pts[5],pts[6])
	lines[4] = geom.add_line(pts[6],pts[1])
	lines[5] = geom.add_line(pts[5],pts[3])
	lines[6] = geom.add_line(pts[3],pts[4])
	lines[7] = geom.add_line(pts[4],pts[6])
	# add loop_line
	line_loops = {}
	line_loops[1] = geom.add_line_loop([lines[1],lines[2],lines[3],lines[4]])
	line_loops[2] = geom.add_line_loop([-lines[3],lines[5],lines[6],lines[7]])
	# add surface
	surfs = {}
	surfs[1] = geom.add_plane_surface(line_loops[1])
	surfs[2] = geom.add_plane_surface(line_loops[2])
	# add pysical_line
	physicalLines = {}
	physicalLines[1] = geom.add_physical_line([lines[1]],label='bd11')
	physicalLines[2] = geom.add_physical_line([lines[6]],label='bd12')
	physicalLines[3] = geom.add_physical_line([lines[2],lines[5]],label='bd21')
	physicalLines[4] = geom.add_physical_line([lines[4],lines[7]],label='bd31')
	# add pysical_surface
	phySurfs = {}
	phySurfs[1] = geom.add_physical_surface([surfs[1]],label='mat1')
	phySurfs[2] = geom.add_physical_surface([surfs[2]],label='mat2')
	points, cells, point_data, cell_data, field_data = pg.generate_mesh(geom)
	# meshing
	mesh = meshio.Mesh(points,cells,point_data,cell_data,field_data)
	# write to disk
	meshio.write(meshPath, mesh)
	# print(geom.get_code())
	return mesh
	
def getBoundaries(mesh,bdParams):
	import numpy as np
	bds = {}
	# lines in each class of boundary
	for keyDict in bdParams.keys():
		if keyDict == 'bd1':
			alpha = 0.
		elif keyDict == 'bd2':
			alpha = 0.
		tempDict = {}
		for key in bdParams[keyDict].keys():
			keyNodes = np.array(mesh.cells['line'][mesh.cell_data['line']['gmsh:physical']==mesh.field_data[key][0]])
			keyNodes = np.unique(keyNodes)			
			tempDict[key+'Node'] = keyNodes
			if keyDict == 'bd3':
				alpha = bdParams[keyDict][key][0]
				tempDict[key+'Params'] = [alpha,alpha*bdParams[keyDict][key][1]]
			else:
				tempDict[key+'Params'] = [alpha,bdParams[keyDict][key]]
		bds[keyDict] = tempDict
	# return a dictionary as the follows 
	# return {'bd1': {'bdNode11': bdNode11, 'bdT11': [0., T1], 
					# 'bdNode12': bdNode12, 'bdT12': [0., T2]}, 
			# 'bd2': {'bdNode21': bdNode21, 'bdq21': [0., q]}, 
			# 'bd3': {'bdEx31': [alpha, alpha*Ts]}   # [alpha, alpha*Ts]
			# }
	# bd1--1st, T0 
	# bd2--2nd, q and adiabatic (q=0)
	# bd3--3rd, alpha, Ts 
	return bds
	
	
def getNodes(mesh):
	return(mesh.points)

def getSurfsWithTrisAndParams(mesh,matParams):
	import numpy as np
	# {'mat1':{tris:[],params:[kappa,miu,miuW,vx,vy,Q]},'mat2':{...}}
	mats = {}
	for key in matParams.keys():
		mats[key] = {'tris':mesh.cells['triangle'][mesh.cell_data['triangle']['gmsh:physical']==mesh.field_data[key][0]],
					'params':matParams[key]}
	return(mats)