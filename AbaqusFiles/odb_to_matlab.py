"""
odb_to_matlab.py
----------------
Extracts modal analysis results from an Abaqus ODB file and saves them
to a MATLAB .mat file containing:
  - phi   : mode shape matrix  [n_dof x n_modes]  sorted by node label & DOF
  - fn    : natural frequencies [n_modes x 1]  (Hz)
  - dof   : DOF table  [n_dof x 1]  (node_label.dof_index 1-6)
  - nodes : node coordinates  [n_nodes x 4]  (label | X | Y | Z)
  - elems : element connectivity  (label | VTK | n1 | n2 | ...)

Usage (must be run inside Abaqus Python):
  abaqus python odb_to_matlab.py --odb path/to/result.odb [options]

Options:
  --step     Step name            (default: last step)
  --instance Assembly instance    (default: first instance)
  --mat      Output .mat path     (default: <odb_name>.mat)

Dependencies:
  odbAccess  - bundled with Abaqus
  numpy      - bundled with Abaqus Python
  scipy      - pip install scipy  (for savemat)
"""

import sys
import os
import argparse
import numpy as np

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Export Abaqus modal ODB to MATLAB .mat")
parser.add_argument("--odb",      required=True, help="Path to the .odb file")
parser.add_argument("--step",     default=None,  help="Step name (default: last step)")
parser.add_argument("--mat",      default=None,  help="Output .mat file")
parser.add_argument("--instance", default=None,  help="Assembly instance name")
args = parser.parse_args()

odb_path = os.path.abspath(args.odb)
mat_path = args.mat or os.path.splitext(odb_path)[0] + ".mat"

# ---------------------------------------------------------------------------
# Open ODB
# ---------------------------------------------------------------------------
try:
    from odbAccess import openOdb
except ImportError:
    sys.exit(
        "ERROR: odbAccess not found.\n"
        "Run this script with:  abaqus python odb_to_matlab.py --odb <file.odb>"
    )

print("Opening ODB: {}".format(odb_path))
odb = openOdb(path=odb_path, readOnly=True)

# ---------------------------------------------------------------------------
# Select step
# ---------------------------------------------------------------------------
step_names = list(odb.steps.keys())
if args.step:
    if args.step not in odb.steps:
        odb.close()
        sys.exit("ERROR: Step '{}' not found. Available: {}".format(args.step, step_names))
    step = odb.steps[args.step]
else:
    step = odb.steps[step_names[-1]]
    print("Using step: '{}' (last step)".format(step.name))

# ---------------------------------------------------------------------------
# Select instance
# ---------------------------------------------------------------------------
instance_names = list(odb.rootAssembly.instances.keys())
if args.instance:
    if args.instance not in odb.rootAssembly.instances:
        odb.close()
        sys.exit("ERROR: Instance '{}' not found. Available: {}".format(
            args.instance, instance_names))
    instance = odb.rootAssembly.instances[args.instance]
else:
    instance = odb.rootAssembly.instances[instance_names[0]]
    print("Using instance: '{}'".format(instance.name))

# ---------------------------------------------------------------------------
# Extract nodes  ->  [n_nodes x 4]  (label, X, Y, Z)
# ---------------------------------------------------------------------------
print("Extracting nodes ...")
node_labels = []
node_coords = []
for node in instance.nodes:
    node_labels.append(node.label)
    coords = list(node.coordinates)
    while len(coords) < 3:       # pad to 3D for 2D models
        coords.append(0.0)
    node_coords.append(coords)

node_labels = np.array(node_labels, dtype=np.int32)
node_coords = np.array(node_coords, dtype=np.float64)

# Sort by node label so everything is in a predictable order
sort_idx    = np.argsort(node_labels)
node_labels = node_labels[sort_idx]
node_coords = node_coords[sort_idx]
nodes_mat   = np.column_stack([node_labels, node_coords])   # [n_nodes x 4]

# Map node label -> sorted row index (0-based)
label_to_idx = {int(lbl): i for i, lbl in enumerate(node_labels)}
n_nodes = len(node_labels)
print("  {} nodes found.".format(n_nodes))

# ---------------------------------------------------------------------------
# Extract elements  ->  single [n_elems x 10] matrix
#   col 0      : element label (ID)
#   col 1      : VTK element type code
#   cols 2-9   : node connectivity (zero-padded for elements with < 8 nodes)
#
# VTK type code reference (most common Abaqus types):
#   Solids : C3D4=10  C3D6=13  C3D8/R=12  C3D10=24  C3D15=26  C3D20/R=25
#   Shells : S3/R=5   S4/R=9   S8R=23
#   Beams  : B31=3    B32=21
#   Points : MASS=1   ROTARYI=1
# ---------------------------------------------------------------------------
VTK_TYPE = {
    # --- Continuum solids ---
    "C3D4"   : 10, "C3D4H"  : 10,
    "C3D6"   : 13, "C3D6H"  : 13,
    "C3D8"   : 12, "C3D8R"  : 12, "C3D8H"  : 12, "C3D8RH" : 12,
    "C3D10"  : 24, "C3D10H" : 24, "C3D10M" : 24, "C3D10MH": 24,
    "C3D15"  : 26, "C3D15H" : 26,
    "C3D20"  : 25, "C3D20R" : 25, "C3D20H" : 25, "C3D20RH": 25,
    # --- Shells ---
    "S3"     :  5, "S3R"    :  5, "STRI3"  :  5,
    "S4"     :  9, "S4R"    :  9, "S4R5"   :  9,
    "S8R"    : 23, "S8R5"   : 23,
    # --- Beams ---
    "B31"    :  3, "B31H"   :  3,
    "B32"    : 21, "B32H"   : 21,
    # --- Points ---
    "MASS"   :  1, "ROTARYI":  1,
}
N_CONNECTIVITY = 8   # fixed connectivity columns (cols 2-9)
N_COLS         = 10  # 1 (label) + 1 (vtk) + 8 (connectivity)

print("Extracting elements ...")
elem_rows     = []
unknown_types = set()
type_counts   = {}

for elem in instance.elements:
    etype = elem.type
    vtk   = VTK_TYPE.get(etype, -1)
    if vtk == -1:
        unknown_types.add(etype)

    conn = list(elem.connectivity)
    conn = conn[:N_CONNECTIVITY] + [0] * max(0, N_CONNECTIVITY - len(conn))
    elem_rows.append([elem.label, vtk] + conn)
    type_counts[etype] = type_counts.get(etype, 0) + 1

elems_mat = np.array(elem_rows, dtype=np.int32)                     # [n_elems x 10]
# Sort: primary = VTK type descending, secondary = element label ascending
sort_idx  = np.lexsort((elems_mat[:, 0], -elems_mat[:, 1]))
elems_mat = elems_mat[sort_idx]

total_elems = len(elems_mat)
for etype, count in sorted(type_counts.items()):
    vtk = VTK_TYPE.get(etype, -1)
    print("  {:>7} elements  type={:<10}  VTK={}".format(count, etype, vtk))
print("  {:>7} elements total.".format(total_elems))

if unknown_types:
    print("  WARNING: Unknown element types (VTK code set to -1): {}".format(
        ", ".join(sorted(unknown_types))))
    print("  Add them to the VTK_TYPE dict at the top of the script.")

# ---------------------------------------------------------------------------
# Extract modal frequencies and mode shapes (bulk read)
# ---------------------------------------------------------------------------
print("Extracting frequencies and mode shapes ...")

# DOF layout per node: U1 U2 U3 UR1 UR2 UR3  (Abaqus indices 1-6)
N_DOF_PER_NODE = 6
U_DOFS         = [0, 1, 2]   # column index in U  data -> DOF slots 0,1,2
UR_DOFS        = [3, 4, 5]   # column index in UR data -> DOF slots 3,4,5

frames      = step.frames
mode_frames = [f for f in frames if f.frameValue > 0]   # skip frame 0 (base state)
n_modes     = len(mode_frames)
print("  {} modes found.".format(n_modes))

# Pre-allocate full phi matrix  [n_nodes*6 x n_modes]
phi = np.zeros((n_nodes * N_DOF_PER_NODE, n_modes), dtype=np.float64)
fn  = np.zeros(n_modes, dtype=np.float64)

# Build DOF table aligned to the sorted node order.
# Rows are ordered: node0/DOF1 ... node0/DOF6, node1/DOF1 ... nodeN/DOF6
# This is kept consistent with how phi rows are written during bulk extraction.
dof_node  = np.repeat(node_labels, N_DOF_PER_NODE)
dof_index = np.tile(np.arange(1, N_DOF_PER_NODE + 1, dtype=np.int32), n_nodes)
# Full 2-col table used internally to track phi row layout
dof_full  = np.column_stack([dof_node, dof_index])      # [n_dof x 2]

for m_idx, frame in enumerate(mode_frames):

    # --- Parse frequency from frame description --------------------------------
    # Abaqus description looks like:
    #   "Mode   1: Value =  1.23456E+04  Freq =  1.7684E+01  (cycles/time)"
    # frame.frameValue is just the mode number (1, 2, 3...), NOT the frequency.
    freq_hz = None
    desc    = frame.description
    if "Freq" in desc:
        try:
            freq_hz = float(desc.split("Freq")[1].split()[1])
        except Exception:
            freq_hz = None
    if freq_hz is None and "Value" in desc:
        try:
            eigenvalue = float(desc.split("Value")[1].split()[1])
            freq_hz    = np.sqrt(abs(eigenvalue)) / (2.0 * np.pi)
        except Exception:
            freq_hz = None
    if freq_hz is None:
        freq_hz = 0.0
        print("    WARNING: Could not parse frequency for mode {}, set to 0.".format(m_idx + 1))
    fn[m_idx] = freq_hz
    print("  Mode {:>3}: {:>12.4f} Hz  | {}".format(m_idx + 1, fn[m_idx], desc.strip()))

    # --- Translational DOFs (U) via bulk read ----------------------------------
    if "U" not in frame.fieldOutputs:
        print("    WARNING: 'U' field not found in mode {}, skipping.".format(m_idx + 1))
        continue

    for block in frame.fieldOutputs["U"].bulkDataBlocks:
        # block.nodeLabels : int array  [n_block]
        # block.dataDouble : float64    [n_block x n_components]
        try:
            data = block.dataDouble
        except Exception:
            data = np.array(block.data, dtype=np.float64)

        labels     = block.nodeLabels
        valid_mask = np.array([int(lbl) in label_to_idx for lbl in labels])
        if not np.any(valid_mask):
            continue

        row_indices = np.array([label_to_idx[int(lbl)]
                                 for lbl in labels[valid_mask]], dtype=np.int32)
        data = data[valid_mask]

        for comp, dof_slot in enumerate(U_DOFS):
            if comp >= data.shape[1]:
                break
            phi[row_indices * N_DOF_PER_NODE + dof_slot, m_idx] = data[:, comp]

    # --- Rotational DOFs (UR) via bulk read ------------------------------------
    if "UR" in frame.fieldOutputs:
        for block in frame.fieldOutputs["UR"].bulkDataBlocks:
            try:
                data = block.dataDouble
            except Exception:
                data = np.array(block.data, dtype=np.float64)

            labels     = block.nodeLabels
            valid_mask = np.array([int(lbl) in label_to_idx for lbl in labels])
            if not np.any(valid_mask):
                continue

            row_indices = np.array([label_to_idx[int(lbl)]
                                     for lbl in labels[valid_mask]], dtype=np.int32)
            data = data[valid_mask]

            for comp, dof_slot in enumerate(UR_DOFS):
                if comp >= data.shape[1]:
                    break
                phi[row_indices * N_DOF_PER_NODE + dof_slot, m_idx] = data[:, comp]

odb.close()
print("ODB closed.")

# ---------------------------------------------------------------------------
# phi: keep full 6-DOF zero-padded layout  [n_nodes*6 x n_modes]
# dof: active DOFs only, single column, format node.dof  (e.g. 1042.2)
# ---------------------------------------------------------------------------
# phi retains a row for every DOF (active or not). Inactive DOFs (e.g.
# rotational DOFs on solid nodes) are zero rows.
phi_reduced = phi

# Active mask: rows that are non-zero in at least one mode
active_mask = np.any(phi != 0, axis=1)

# Build dof vector: format node_label.dof_index (e.g. node 1042, DOF 2 -> 1042.2)
# Use float64 so the decimal digit is preserved exactly (dof index is 1-6,
# so the fractional part is always a single digit with no rounding issues).
dof_active = (dof_full[active_mask, 0].astype(np.float64)
              + dof_full[active_mask, 1].astype(np.float64) * 0.1)
dof_active = dof_active.reshape(-1, 1)   # column vector for MATLAB

print("phi shape : {} (6 DOFs per node, zeros for inactive DOFs).".format(phi.shape))
print("dof shape : {} (active DOFs only, format = node.dof).".format(dof_active.shape))
print("DOF ordering: sorted by ascending node label, then DOF index 1->6.")

# ---------------------------------------------------------------------------
# Save to .mat
# ---------------------------------------------------------------------------
try:
    from scipy.io import savemat
except ImportError:
    sys.exit(
        "ERROR: scipy not available.\n"
        "Install with:  pip install scipy"
    )

save_dict = {
    "phi"  : phi_reduced,          # [n_nodes*6 x n_modes]
    "fn"   : fn.reshape(-1, 1),    # [n_modes x 1]  column vector in MATLAB
    "dof"  : dof_active,           # [n_active_dof x 1]  format: node.dof
    "nodes": nodes_mat,            # [n_nodes x 4]  (label | X | Y | Z)
    "elems": elems_mat,            # [n_elems x 10] (label | vtk_type | n1..n8)
}

savemat(mat_path, save_dict, do_compression=True)

print("\nSaved: {}".format(mat_path))
print("\nVariables in .mat file:")
print("  phi   {} - mode shape matrix (6 DOFs/node, zero-padded)".format(phi_reduced.shape))
print("  fn    {} - natural frequencies (Hz)".format(fn.reshape(-1, 1).shape))
print("  dof   {} - active DOFs, format node.dof (e.g. 1042.2 = node 1042, U2)".format(dof_active.shape))
print("  nodes {} - (label | X | Y | Z)".format(nodes_mat.shape))
print("  elems {} - (label | vtk_type | n1..n8)".format(elems_mat.shape))
print("\nElement matrix column layout:")
print("  col 1    : element label")
print("  col 2    : VTK type code  (1=point 3=beam 5=tri 9=quad 10=tet 12=hex ...)")
print("  col 3-10 : node connectivity (zero-padded for elements with < 8 nodes)")
print("\nDOF index legend:  1=U1  2=U2  3=U3  4=UR1  5=UR2  6=UR3")
print("\nDone.")