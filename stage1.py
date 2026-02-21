# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 12:23:07 2025

@author: dkummara-yal
"""

# -----------------------------------------------
# Abaqus Script
# -----------------------------------------------
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from abaqusConstants import *
from abaqus import session
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import time
import math
import numpy as np
from odbAccess import openOdb


# Parameters according to the research paper

# Parameters according to the research paper

ht =60/2
r_coil=29
r_out=32/2
nr = 3             
ns=30
nc=16
r_in=23/2


# my Parameters

outer_diameter = r_out*2            # mm, Outer diameter of wire
inner_diameter = r_in*2            # mm, Inner diameter of wire
mean_diameter  = r_coil *2           # mm, Mean diameter of the coil
pitch_height   = ht*2           # mm, Pitch of the coil
number_of_turns = 0.5          # Number of turns

mesh_seed_main = int(nc/2)             # Mesh seed for main edges
mesh_seed_medium = ns             # Mesh seed for medium edges


no_of_sections = ns
coil_radius  = r_coil
outer_radius = r_out
inner_radius = r_in
revolve_angle = number_of_turns * 360.0
y0 = 0           # Hauteur de coupe de fibre neutre en 0
yF = ht    # Hauteur de coupe de fibre neutre en F
U = 46 



start_time = time.time()

# Naming
ModelName = 'Model-1'
PartName = 'Part-1'
InstanceName = PartName + '-1'
jobName = 'step1_spring'


# Create model and part
myModel = mdb.Model(name=ModelName)

# Sketch hollow circular wire profile
s = myModel.ConstrainedSketch(name='__profile__', sheetSize=200.0)
s.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
s.FixedConstraint(entity=s.geometry[2])
s.ConstructionLine(point1=(0.0, 30.0), point2=(0.0, -30.0))
s.VerticalConstraint(addUndoState=False, entity=s.geometry[3])
# Outer circle (wire outer diameter)
s.CircleByCenterPerimeter(center=(coil_radius, 0.0), point1=(coil_radius + outer_radius, 0.0))
# Inner circle (wire inner diameter - hollow)
s.CircleByCenterPerimeter(center=(coil_radius, 0.0), point1=(coil_radius + inner_radius, 0.0))
s.sketchOptions.setValues(constructionGeometry=ON)
s.assignCenterline(line=s.geometry[2])

# Create part and revolve to make 3D solid spring
myPart = myModel.Part(name=PartName, dimensionality=THREE_D, type=DEFORMABLE_BODY)
myPart.BaseSolidRevolve(
    angle=revolve_angle,
    flipPitchDirection=OFF,
    flipRevolveDirection=OFF,
    moveSketchNormalToPath=OFF,
    pitch=pitch_height,
    sketch=s
)

# Delete sketch to clean up
del myModel.sketches['__profile__']

# Partition faces as per original code
myPart.PartitionFaceByShortestPath(
    faces=myPart.faces.getSequenceFromMask(('[#4 ]', ), ),
    point1=myPart.vertices[2],
    point2=myPart.vertices[0]
)

myPart.PartitionFaceByShortestPath(
    faces=myPart.faces.getSequenceFromMask(('[#4 ]', ), ),
    point1=myPart.InterestingPoint(myPart.edges[0], MIDDLE),
    point2=myPart.InterestingPoint(myPart.edges[3], MIDDLE)
)

    
# --- ASSEMBLY

# Creation de l assemblage
myAssembly = myModel.rootAssembly
myAssembly.DatumCsysByDefault(CARTESIAN)

# Ajout d une instance dans l assemblage
myInstance = myAssembly.Instance(dependent=ON, name=InstanceName, part=myPart)

# Creation des plans en y0 et yF correspondants aux cylindres supports
plane0 = myPart.DatumPlaneByPrincipalPlane(offset=y0, principalPlane=XZPLANE)
planeF = myPart.DatumPlaneByPrincipalPlane(offset=-yF, principalPlane=XZPLANE)

# --- STEP

# Creation du step
myModel.StaticStep(name='Step-1', previous='Initial') # Rajouter un ON pour non linearite geometry

# Variables en sortie du Job
# S contrainte sur les 4 points d integration, U deplacement, MISESMAX contrainte max en tout point d integration, RF effort associe
myModel.FieldOutputRequest(name='F-Output-1', createStepName='Step-1', variables=('S','U','MISESMAX','RF'))



# Seeding and meshing
myPart.seedEdgeByNumber(constraint=FINER, edges=myPart.edges.getSequenceFromMask(('[#82 ]', ), ), number=mesh_seed_main)
myPart.seedEdgeByNumber(constraint=FINER, edges=myPart.edges.getSequenceFromMask(('[#18 ]', ), ), number=mesh_seed_main)
myPart.seedEdgeByNumber(constraint=FINER, edges=myPart.edges.getSequenceFromMask(('[#5 ]', ), ), number=nr)
myPart.seedEdgeByNumber(constraint=FINER, edges=myPart.edges.getSequenceFromMask(('[#444420 ]', ), ), number=mesh_seed_medium)
myPart.seedEdgeByNumber(constraint=FINER, edges=myPart.edges.getSequenceFromMask(('[#11100 ]', ), ), number=mesh_seed_medium)
myPart.generateMesh()


myJob = mdb.Job(name=jobName, model=ModelName) 

myJob.writeInput()

