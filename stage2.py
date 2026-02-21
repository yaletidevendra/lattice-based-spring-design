# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 12:21:23 2025

@author: dkummara-yal
"""


import re
import numpy as np
import pyxel as px
import matplotlib.pyplot as plt
from collections import defaultdict

plt.close('all')


# -------------------------------
# Structure class
# -------------------------------
class Structure:
    def __init__(self, name, catalog_content):
        self.name = name.lower()
        self.catalog_content = catalog_content
        self.nodes = None
        self.elements = None
        self._load_structure_from_catalog()

    def _load_structure_from_catalog(self):
        blocks = self.catalog_content.split(
            "-----------------------------------------------------------------------------------------"
        )

        for block in blocks:
            name_line_match = re.search(r"Name:\s*(.+)", block, re.IGNORECASE)
            if not name_line_match:
                continue
            found_name = name_line_match.group(1).strip().lower()
            if self.name == found_name:
                nodes_match = re.search(r"Nodal positions:\s*\n([\d\.\-\s]+)", block, re.IGNORECASE)
                elems_match = re.search(r"Bar connectivities:\s*\n([\d\s\n]+)", block, re.IGNORECASE)

                if not nodes_match or not elems_match:
                    raise ValueError(f"Could not find nodes or elements for '{self.name}'")

                node_lines = nodes_match.group(1).strip().splitlines()
                self.nodes = np.array([list(map(float, line.strip().split())) for line in node_lines]) * 2 - 1

                elem_lines = elems_match.group(1).strip().splitlines()
                self.elements = {1: np.array([list(map(lambda x: int(x)-1, line.strip().split())) for line in elem_lines])}
                return

        raise ValueError(f"Structure '{self.name}' not found in catalog")

    def rotated(self, axis=None, angle_deg=90):
        if axis is None:
            return self.nodes.copy()

        theta = np.radians(angle_deg)
        if axis == 'y':
            R = np.array([
                [np.cos(theta), 0, np.sin(theta)],
                [0, 1, 0],
                [-np.sin(theta), 0, np.cos(theta)]
            ])
        elif axis == 'z':
            R = np.array([
                [np.cos(theta), -np.sin(theta), 0],
                [np.sin(theta),  np.cos(theta), 0],
                [0, 0, 1]
            ])
        elif axis == 'x':
            R = np.array([
                [1, 0, 0],
                [0, np.cos(theta), -np.sin(theta)],
                [0, np.sin(theta),  np.cos(theta)]
            ])
        else:
            raise ValueError(f"Invalid axis '{axis}'. Choose from 'x', 'y', or 'z'.")
        return self.nodes @ R.T


# -------------------------------
# MeshComposer class
# -------------------------------
class MeshComposer:
    def __init__(self, base_mesh_file):
        self.base_mesh = px.ReadMesh(base_mesh_file, 3)
        self.output_files = []
        self.r_min_values = []

    @staticmethod
    def compute_element_lengths(mg):
        nodes = mg.n
        element_type = list(mg.e.keys())[0]
        elements = mg.e[element_type]
        lengths = np.linalg.norm(nodes[elements[:, 1]] - nodes[elements[:, 0]], axis=1)
        return lengths

    def compose_with_structure(self, structure: Structure, rotation_cases=None):
        if rotation_cases is None:
            rotation_cases = {"none": None}

        for suffix, axis in rotation_cases.items():
            rotated_nodes = structure.rotated(axis=axis)
            new_mesh = px.Mesh(structure.elements, rotated_nodes, dim=3)
            # new_mesh.Plot()

            composed_mesh = self.base_mesh.FEComposition(new_mesh)

            lengths = self.compute_element_lengths(composed_mesh)
            L_min, L_avg, L_max = lengths.min(), lengths.mean(), lengths.max()
            r_min, r_avg, r_max = L_min / 16, L_avg / 16, L_max / 16

            print(f"\nStructure: {structure.name}, Rotation: {suffix}")
            print(f"  L_min = {L_min:.6f}, L_avg = {L_avg:.6f}, L_max = {L_max:.6f}")
            print(f"  r_min = {r_min:.6f}, r_avg = {r_avg:.6f}, r_max = {r_max:.6f}")

            outname = f"{structure.name.replace('.', '_')}-{suffix}.inp"
            composed_mesh.Write(outname)
            self.output_files.append(outname)
            self.r_min_values.append(r_min)
            print(f"Wrote: {outname}")

    def save_output_list(self, filename):
        with open(filename, "w") as f:
            for name, rmin in zip(self.output_files, self.r_min_values):
                f.write(f"{name}    {rmin:.6f}\n")
        print(f"All compositions completed and saved to '{filename}'")


# -------------------------------
# Main script
# -------------------------------
if __name__ == "__main__":
    # Load catalog
    with open("Unit_Cell_Catalog.txt", "r") as f:
        catalog_content = f.read()

    # Group by normalized unit cell parameters
    group_dict = defaultdict(list)
    blocks = catalog_content.split(
        "-----------------------------------------------------------------------------------------"
    )

    for block in blocks:
        name_match = re.search(r"Name:\s*(.+)", block, re.IGNORECASE)
        params_match = re.search(
            r"Normalized unit cell parameters\s*\(a,b,c,alpha,beta,gamma\):\s*\n([\d\.,\s]+)",
            block,
            re.IGNORECASE
        )
        if not name_match or not params_match:
            continue

        name = name_match.group(1).strip()
        params_line = params_match.group(1).strip()
        params = tuple(map(float, params_line.split(',')))
        group_dict[params].append(name)

    # # Show all groups and counts
    # print("Unit cell groups and counts:")
    # for params, names in group_dict.items():
    #     print(f"Parameters: {params}, Count: {len(names)}")

    # --- Select a specific group ---
    target_params = (1.0, 1.0, 1.0, 90.0, 90.0, 90.0)  # example cubic
    if target_params not in group_dict:
        print("No unit cells found with these parameters.")
        exit()

    selected_cells = group_dict[target_params]
    print(f"Total cells in selected group: {len(selected_cells)}")

    # --- Select range to process ---
    start_idx, end_idx = 0, 325  # first 10 cells
    cells_to_generate = selected_cells[start_idx:end_idx]

    # --- Generate INP files and compute r_min ---
    composer = MeshComposer('step1_spring.inp')
    rotation_cases = {"none": None}

    for name in cells_to_generate:
        print(f"\nProcessing structure: {name}")
        structure = Structure(name, catalog_content)
        composer.compose_with_structure(structure, rotation_cases)

    # --- Save results ---
    composer.save_output_list("inp_selected_rmin.txt")