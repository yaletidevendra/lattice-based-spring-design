# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 12:21:56 2025

@author: dkummara-yal
"""



print("\n" * 5)  # Pushes old output up
print("========== Batch Run Start ==========")

from abaqus import *
from abaqusConstants import *
import time
import numpy as np
import os
import re
from odbAccess import openOdb

# Read input file paths and corresponding r_beam values from the text file
input_list_file = r"C:/ABAQUS_SCRATCH/a1/inp_selected_rmin.txt"  # <-- Your txt file here

input_files = []
r_beams = []

with open(input_list_file, "r") as f:
    for line in f:
        if line.strip():
            parts = line.strip().split()
            input_files.append(parts[0])       # file name
            r_beams.append(float(parts[1]))    # r_beam

# Write header only once
results_file = "C:/ABAQUS_SCRATCH/a1/results_summary.txt"
if not os.path.exists(results_file):
    with open(results_file, "w") as result_file:
        result_file.write("jobName   r_beam(mm)   S(MPa)   F(N)   U(mm)   K(N/mm)   K/mass(N/mm/kg)   mass(kg)\n")
        
# Loop over all input files and run simulations
for idx, (inputFileName, r_beam) in enumerate(zip(input_files, r_beams), start=1):
    print("\n" * 5)
    print("========== Starting Run {}/{} ==========".format(idx, len(input_files)))
    print("Running simulation on input file: {}".format(inputFileName))

    start_time = time.time()

    try:
        # -------------------------------
        # PARAMETERS
        # -------------------------------
        ht = 60/2
        r_coil = 29
        r_out = 32/2
        nr = 3
        ns = 30
        nc = 16
        r_in = 23/2

        # Derived Parameters
        outer_diameter = r_out*2
        inner_diameter = r_in*2
        mean_diameter  = r_coil *2
        pitch_height   = ht*2
        number_of_turns = 0.5

        mesh_seed_main = int(nc/2)
        mesh_seed_medium = ns

        no_of_sections = ns
        coil_radius  = r_coil
        outer_radius = r_out
        inner_radius = r_in
        revolve_angle = number_of_turns * 360.0
        y0 = -ht
        yF = 0

        mesh_seed_small = nr
        MeshSize = 2.25
        F = 0.5  # Applied force (N)

        # -------------------------------
        # UNIQUE NAMES
        # -------------------------------
        filename = os.path.basename(inputFileName)
        base_name = filename.replace('.inp', '')
        base_name = re.sub(r'\W+', '_', base_name)
        if base_name[0].isdigit():
            base_name = '_' + base_name
        ModelName = 'Model_{}'.format(base_name)
        PartName = 'PART-1'
        jobName = '{}'.format(base_name)
        InstanceName = PartName + '-1'

        # -------------------------------
        # MODEL IMPORT
        # -------------------------------
        mdb.ModelFromInputFile(inputFileName=inputFileName, name=ModelName)
        myModel = mdb.models[ModelName]
        myPart = myModel.parts[PartName]

        # Remove existing instance if it exists
        if 'PART-1-1' in myModel.rootAssembly.features.keys():
            del mdb.models[ModelName].rootAssembly.features['PART-1-1']

        # -------------------------------
        # MATERIAL + SECTION
        # -------------------------------
        myMaterial = myModel.Material(name='Nylon-PA12')
        myMaterial.Elastic(table=((1580, 0.38), ))
        myMaterial.Density(table=((1.033e-9, ), ))

        myProfile = myModel.CircularProfile(name='Profile-1', r=r_beam)  # <-- use r_beam from txt
        myModel.BeamSection(consistentMassMatrix=False, integration=DURING_ANALYSIS,
                            material='Nylon-PA12', name='Section-1', poissonRatio=0.0,
                            profile='Profile-1', temperatureVar=LINEAR)

        ALL = myPart.Set(name='ALL', elements=myPart.elements)
        myPart.assignBeamSectionOrientation(method=N1_COSINES, n1=(1.0, 1.0, 0.0), region=ALL)
        myPart.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE,
                                 region=ALL, sectionName='Section-1', thicknessAssignment=FROM_SECTION)

        # -------------------------------
        # ASSEMBLY
        # -------------------------------
        myAssembly = myModel.rootAssembly
        myAssembly.DatumCsysByDefault(CARTESIAN)
        myInstance = myAssembly.Instance(dependent=ON, name=InstanceName, part=myPart)

        plane0 = myPart.DatumPlaneByPrincipalPlane(offset=y0, principalPlane=XZPLANE)
        planeF = myPart.DatumPlaneByPrincipalPlane(offset=yF, principalPlane=XZPLANE)

        myAssembly.regenerate()
        mass_props = myAssembly.getMassProperties(regions=(myInstance,))
        mass = mass_props['mass']

        # -------------------------------
        # STEP
        # -------------------------------
        myModel.StaticStep(name='Step-1', previous='Initial', initialInc=0.01,
                           minInc=1e-6, maxNumInc=10000)
        myModel.FieldOutputRequest(name='F-Output-1', createStepName='Step-1',
                                   variables=('S', 'U', 'MISESMAX', 'RF'))
        myModel.HistoryOutputRequest(createStepName='Step-1', name='H-Output-2',
                                     variables=('MASS', ))

        # -------------------------------
        # INTERACTIONS + LOAD
        # -------------------------------
        pt0 = myAssembly.ReferencePoint(point=(0.0, y0, 0.0))
        ptF = myAssembly.ReferencePoint(point=(0.0, yF, 0.0))
        RPF = myAssembly.Set(name='m_Set-f', referencePoints=(myAssembly.referencePoints[ptF.id], ))
        RP0 = myAssembly.Set(name='m_Set-0', referencePoints=(myAssembly.referencePoints[pt0.id], ))

        tolz = 1.0
        nodes_plane_F = myInstance.nodes.getByBoundingBox(
            xMin=coil_radius - outer_radius, xMax=coil_radius + outer_radius,
            yMin=yF - outer_radius, yMax=yF + outer_radius,
            zMin=-tolz, zMax=+tolz
        )
        myAssembly.Set(name='AllVerticesOnPlaneF', nodes=nodes_plane_F)

        nodes_plane_0 = myInstance.nodes.getByBoundingBox(
            xMin=-coil_radius - outer_radius, xMax=-coil_radius + outer_radius,
            yMin=y0 - outer_radius, yMax=y0 + outer_radius,
            zMin=-tolz, zMax=tolz
        )
        myAssembly.Set(name='AllVerticesOnPlane0', nodes=nodes_plane_0)

        myModel.Coupling(controlPoint=RP0, couplingType=KINEMATIC, influenceRadius=WHOLE_SURFACE,
                         localCsys=None, name='Constraint-0', surface=myModel.rootAssembly.sets['AllVerticesOnPlane0'],
                         u1=ON, u2=ON, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)

        myModel.Coupling(controlPoint=RPF, couplingType=KINEMATIC, influenceRadius=WHOLE_SURFACE,
                         localCsys=None, name='Constraint-F', surface=myModel.rootAssembly.sets['AllVerticesOnPlaneF'],
                         u1=OFF, u2=ON, u3=OFF, ur1=OFF, ur2=OFF, ur3=OFF)

        myModel.EncastreBC(createStepName='Initial', name='BC-0', region=RP0)
        myModel.ConcentratedForce(name='Load-F', createStepName='Step-1', region=RPF, cf2=-F,
                                  distributionType=UNIFORM, field='')

        # -------------------------------
        # JOB
        # -------------------------------
        myJob = mdb.Job(name=jobName, model=ModelName)
        myJob.submit(consistencyChecking=OFF)
        myJob.waitForCompletion()

        # -------------------------------
        # RESULTS
        # -------------------------------
        odb = openOdb(path=jobName+'.odb')
        frame = odb.steps['Step-1'].frames[-1]
        region = odb.steps['Step-1'].historyRegions['Assembly ASSEMBLY']
        mass_history = region.historyOutputs['MASS']
        mass = mass_history.data[-1][1]
        mass_kg = mass * 1000

        misesmax_field = frame.fieldOutputs['MISESMAX']
        max_S = max(misesmax_field.values, key=lambda misesmax: misesmax.data).data

        disp_field = frame.fieldOutputs['U']
        max_U_value = max(disp_field.values, key=lambda v: (v.data[0]**2 + v.data[1]**2 + v.data[2]**2)**0.5)
        max_U = (max_U_value.data[0]**2 + max_U_value.data[1]**2 + max_U_value.data[2]**2)**0.5

        K = F / max_U
        km = K / mass_kg

        end_time = time.time()
        elapsed_time = end_time - start_time
        print("Finished Run {}/{} in {} seconds.".format(idx, len(input_files), np.round(elapsed_time, 2)))
        print("-----------------------------------------------")
        print('S =', np.round(max_S, 3), 'MPa')
        print('F =', np.round(F, 3), 'N')
        print('U =', np.round(max_U, 3), 'mm')
        print('K =', np.round(K, 5), 'N/mm')
        print("K/mass = {:.6f} N/mm/kg".format(km))
        print("Mass = {:.6e} kg".format(mass_kg))
        print("-----------------------------------------------")
        
        with open(results_file, "a") as result_file:
            result_file.write("{} {:.6f} {:.3f} {:.3f} {:.3f} {:.5f} {:.6f} {:.6e}\n".format(
                jobName, r_beam, max_S, F, max_U, K, km, mass_kg
                ))



        odb.close()

    except Exception as e:
        # If *any* error occurs, log it and continue
        print(" ERROR in Run {}/{}: {}. Skipping...".format(idx, len(input_files), str(e)))
        with open(results_file, "a") as result_file:
            result_file.write("{} ERROR\n".format(os.path.basename(inputFileName)))
        continue

print("\n========== All Runs Completed ==========")