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
	pts[0] = geom.add_point([0,0,0],lcar=1)
	pts[1] = geom.add_point([10,0,0],lcar=1)
	pts[2] = geom.add_point([10,10,0],lcar=1)
	pts[3] = geom.add_point([0,10,0],lcar=1)
	# add line
	lines = {}
	lines[0] = geom.add_line(pts[0],pts[1])
	lines[1] = geom.add_line(pts[1],pts[2])
	lines[2] = geom.add_line(pts[2],pts[3])
	lines[3] = geom.add_line(pts[3],pts[0])
	# add loop_line
	line_loops = {}
	line_loops[0] = geom.add_line_loop([lines[0],lines[1],lines[2],lines[3]])
	# add surface
	surfs = {}
	surfs[0] = geom.add_plane_surface(line_loops[0])
	# add pysical_line
	physicalLines = {}
	physicalLines[0] = geom.add_physical_line([lines[0]],label='bd11')
	physicalLines[1] = geom.add_physical_line([lines[2]],label='bd12')
	physicalLines[2] = geom.add_physical_line([lines[1],lines[3]],label='bd2')
	# add pysical_surface
	phySurfs = {}
	phySurfs[0] = geom.add_physical_surface([surfs[0]],label='phy1')
	points, cells, point_data, cell_data, field_data = pg.generate_mesh(geom)
	# meshing
	mesh = meshio.Mesh(points,cells,point_data,cell_data,field_data)
	# write to disk
	meshio.write(meshPath, mesh)
	# print(geom.get_code())
	return mesh
	
def getBoundaries(mesh,alpha=3.,Ts=0):
	import numpy as np
	# lines in each class of boundary
	bd11 = np.array(mesh.cells['line'][mesh.cell_data['line']['gmsh:physical']==1])
	bd12 = np.array(mesh.cells['line'][mesh.cell_data['line']['gmsh:physical']==2])
	bd21 = np.array(mesh.cells['line'][mesh.cell_data['line']['gmsh:physical']==3])
	# nodes in each class of boundary
	bdNode11 = np.unique(bd11)
	bdNode12 = np.unique(bd12)
	bdNode21 = np.unique(bd21)
	# return a dictionary. 
	# bd1--1st, T0 
	# bd2--2nd, q and adiabatic 
	# bd3--3rd, Ts 
	# [alpha, beta]
	return {'bd1': {'bdNode11': bdNode11, 'bdT11': [1., 10.], 
					'bdNode12': bdNode12, 'bdT12': [1., 100.]}, 
			'bd2': {'bdNode21': bdNode21, 'bdq21': [0., 0.]}, 
			'bd3': {'bdEx31': [alpha, alpha*Ts]}   # [alpha, alpha*Ts]
			}
	
def getNodes(mesh):
	return(mesh.points)

def getTriangles(mesh):
	return(mesh.cells['triangle'])