Overview

This repository contains Python scripts developed during an internship for the automated generation, simulation, and post-processing of lattice-based spring models in Abaqus.

The workflow is divided into three sequential stages, each implemented as a dedicated Python script. Together, these stages enable the creation of lattice-filled spring geometries, execution of finite element simulations, and extraction of mechanical performance metrics.

Workflow Structure
Stage 1: Unit Cell Layout Definition (stage1.py)

This stage defines the layout and discretization of the unit cell within the solid spring geometry.

Stage 2: Parametric Mapping and Lattice Embedding (stage2.py)

This stage performs parametric mapping to replace the solid elements with the selected lattice unit cell.

Stage 3: Abaqus Simulation and Post-Processing (stage3.py)

This stage automates the finite element analysis and result extraction in Abaqus.


Extra: Machine Learning–Based Unit Cell Evaluation (ml.py)

This stage introduces machine learning concepts to efficiently analyze and predict the performance of lattice unit cells.

Additional Information:

Please do not hesitate to contact me if you require:

Full details of the unit cell database

Information on machine learning features and training

Additional scripts or documentation related to this project

📧 Email: yaletidevendra@gmail.com
