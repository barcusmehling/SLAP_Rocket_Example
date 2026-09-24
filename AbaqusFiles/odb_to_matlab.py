"""
odb_to_matlab.py
----------------
Extracts modal analysis results from an Abaqus ODB file and saves them
to a MATLAB .mat file containing:
  - phi      : mode shape matrix     [n_dof x n_modes]   sorted by node label & DOF
  - fn       : natural frequencies   [n_modes x 1]       (Hz)
  - dof      : DOF table             [n_active_dof x 1]  format: node.dof (e.g. 1042.2)
  - nodes    : node coordinates      [n_nodes x 4]       (label | X | Y | Z)
  - elems    : element connectivity  [n_elems x 10]      (label | vtk_type | n1..n8)

  If --elset is provided, also extracts:
  - psi      : modal stress tensor   [n_results x 6 x n_modes]
               component order: S11 S22 S33 S12 S13 S23
  - elset_id : result ID table       [n_results x 2]     (elem_label | int_point)

  If --nset is provided, also extracts:
  - nset     : nodeset vector       [n_nodes x 4]       (label | X | Y | Z)

Usage Examples (Using a Git Bash terminal):
  abaqus python odb_to_matlab.py --odb path/to/result.odb [options]
  abaqus python odb_to_matlab.py --odb path/to/result.odb --elset ELSET_NAME (elset name in all caps)
  abaqus python odb_to_matlab.py --odb path/to/result.odb --nset Nset_Name (nset name case sensitive)
  abaqus python odb_to_matlab.py --odb path/to/result.odb --elset ELSET_NAME --nset Nset_Name
  
Options:
  --step     Step name            (default: last step)
  --instance Assembly instance    (default: first instance)
  --mat      Output .mat path     (default: <odb_name>_Modes.mat)
  --elset    Element set name     (optional, enables stress extraction)
  --nset     Node set name        (optional, enables node set extraction)

Dependencies:
  odbAccess  - bundled with Abaqus Python
  numpy      - bundled with Abaqus Python
  scipy      - bundled with Abaqus Python (not available in some older versions of Abaqus)
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
parser.add_argument("--elset",    default=None,  help="Element set name for stress extraction (optional)")
parser.add_argument("--nset",     default=None,  help="Node set name to be extracted")
args = parser.parse_args()

odb_path    = os.path.abspath(args.odb)
odb_name    = os.path.splitext(os.path.basename(odb_path))[0]
default_dir = os.path.abspath(os.path.join(os.path.dirname(odb_path), "..", "ModeShapes"))
if not os.path.exists(default_dir):
    os.makedirs(default_dir)
mat_path   = args.mat or os.path.join(default_dir, odb_name + "_Modes.mat")
elset_name = args.elset

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
node_coords = np.array(node_coords, dtype=np.float32)

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
# Resolve element set and build stress extraction lookups (if --elset given)
# ---------------------------------------------------------------------------
do_stress    = elset_name is not None
elset_region = None
result_ids   = []
psi          = None

if do_stress:
    print("\nResolving element set '{}' for stress extraction ...".format(elset_name))

    # Check assembly-level sets first, then fall back to per-instance search
    if elset_name in odb.rootAssembly.elementSets:
        elset_region = odb.rootAssembly.elementSets[elset_name]
        print("  Found at assembly level (may span multiple instances).")
    else:
        for iname in instance_names:
            inst = odb.rootAssembly.instances[iname]
            if elset_name in inst.elementSets:
                elset_region = inst.elementSets[elset_name]
                print("  Found on instance '{}'.".format(iname))
                break

    if elset_region is None:
        odb.close()
        sys.exit("ERROR: Element set '{}' not found on any instance or at assembly level.\n"
                 "       Available instances: {}".format(elset_name, instance_names))

    # Build element label -> instance lookup (handles assembly-level elsets)
    elem_to_instance = {}
    elset_labels     = set()
    for el in elset_region.elements:
        elset_labels.add(el.label)
        if hasattr(el, "instanceName"):
            iname = el.instanceName
            elem_to_instance[el.label] = odb.rootAssembly.instances[iname]
        else:
            elem_to_instance[el.label] = instance

    print("  {} elements across {} instance(s).".format(
        len(elset_labels),
        len(set(i.name for i in elem_to_instance.values()))))

    # Deduplicate instances by name
    seen  = set()
    insts = []
    for inst in elem_to_instance.values():
        if inst.name not in seen:
            seen.add(inst.name)
            insts.append(inst)

# ---------------------------------------------------------------------------
# Resolve node set and extract coordinates (if --nset given)
# ---------------------------------------------------------------------------
do_nset   = args.nset is not None
nset_mat  = None

if do_nset:
    nset_name = args.nset
    print("\nResolving node set '{}' ...".format(nset_name))

    nset_region = None
    if nset_name in odb.rootAssembly.nodeSets:
        nset_region = odb.rootAssembly.nodeSets[nset_name]
        print("  Found at assembly level.")
    else:
        for iname in instance_names:
            inst = odb.rootAssembly.instances[iname]
            if nset_name in inst.nodeSets:
                nset_region = inst.nodeSets[nset_name]
                print("  Found on instance '{}'.".format(iname))
                break

    if nset_region is None:
        odb.close()
        sys.exit("ERROR: Node set '{}' not found on any instance or at assembly level.\n"
                 "       Available instances: {}".format(nset_name, instance_names))

    # Extract [node_label, X, Y, Z]  sorted by node label
    nset_rows = []
    for node in nset_region.nodes:
        coords = list(node.coordinates)
        while len(coords) < 3:
            coords.append(0.0)
        nset_rows.append([node.label] + coords)

    nset_rows.sort(key=lambda x: x[0])
    nset_mat = np.array(nset_rows, dtype=np.float32)   # [n_nodes x 4]
    print("  {} nodes found in node set.".format(len(nset_mat)))

# ---------------------------------------------------------------------------
# Extract modal frequencies and mode shapes (bulk read)
# ---------------------------------------------------------------------------
print("\nExtracting frequencies and mode shapes ...")

# DOF layout per node: U1 U2 U3 UR1 UR2 UR3  (Abaqus indices 1-6)
N_DOF_PER_NODE = 6
U_DOFS         = [0, 1, 2]   # column index in U  data -> DOF slots 0,1,2
UR_DOFS        = [3, 4, 5]   # column index in UR data -> DOF slots 3,4,5

frames      = step.frames
mode_frames = [f for f in frames if f.frameValue > 0]   # skip frame 0 (base state)
n_modes     = len(mode_frames)
print("  {} modes found.".format(n_modes))

# Pre-allocate full phi matrix  [n_nodes*6 x n_modes]
phi = np.zeros((n_nodes * N_DOF_PER_NODE, n_modes), dtype=np.float32)
fn  = np.zeros(n_modes, dtype=np.float32)

# Build DOF table aligned to the sorted node order.
# Rows are ordered: node0/DOF1 ... node0/DOF6, node1/DOF1 ... nodeN/DOF6
dof_node  = np.repeat(node_labels, N_DOF_PER_NODE)
dof_index = np.tile(np.arange(1, N_DOF_PER_NODE + 1, dtype=np.int32), n_nodes)
dof_full  = np.column_stack([dof_node, dof_index])      # [n_dof x 2]

# Pre-allocate psi  [n_results x 6 x n_modes]  if stress extraction is on
if do_stress:
    print("\nDiscovering stress result locations from first mode frame ...")
    first_frame  = mode_frames[0]
    subset       = first_frame.fieldOutputs["S"].getSubset(region=elset_region)

    # Collect unique element labels only (average over integration points)
    elem_label_set = []
    seen_els       = set()
    for val in subset.values:
        if val.elementLabel not in seen_els:
            seen_els.add(val.elementLabel)
            elem_label_set.append(val.elementLabel)

    elem_label_set.sort()
    result_ids = np.array(elem_label_set, dtype=np.int32).reshape(-1, 1)  # [n_elems x 1]
    n_results  = len(result_ids)
    id_to_row  = {int(result_ids[i, 0]): i for i in range(n_results)}
    psi        = np.zeros((n_results, 6, n_modes), dtype=np.float32)
    ip_count   = np.zeros(n_results, dtype=np.int32)   # for averaging
    print("  {} elements found in elset.".format(n_results))

# ---------------------------------------------------------------------------
# Main extraction loop — mode shapes and stresses in one pass
# ---------------------------------------------------------------------------
print("\nExtracting mode shapes{} ...".format(
    " and stresses" if do_stress else ""))

for m_idx, frame in enumerate(mode_frames):

    # --- Parse frequency from frame description ----------------------------
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

    # --- Translational DOFs (U) via bulk read ------------------------------
    if "U" not in frame.fieldOutputs:
        print("    WARNING: 'U' field not found in mode {}, skipping.".format(m_idx + 1))
        continue

    for block in frame.fieldOutputs["U"].bulkDataBlocks:
        try:
            data = block.dataDouble
        except Exception:
            data = np.array(block.data, dtype=np.float32)

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

    # --- Rotational DOFs (UR) via bulk read --------------------------------
    if "UR" in frame.fieldOutputs:
        for block in frame.fieldOutputs["UR"].bulkDataBlocks:
            try:
                data = block.dataDouble
            except Exception:
                data = np.array(block.data, dtype=np.float32)

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

    # --- Stresses (S) ------------------------------------------------------
    if do_stress:

        if "S" not in frame.fieldOutputs:
            print("    WARNING: 'S' field not found in mode {}, skipping.".format(m_idx + 1))
        else:
            subset = frame.fieldOutputs["S"].getSubset(region=elset_region)

            # Accumulate stress over integration points
            psi_frame = np.zeros((n_results, 6), dtype=np.float32)
            ip_count  = np.zeros(n_results, dtype=np.int32)

            for val in subset.values:
                el = int(val.elementLabel)
                if el not in id_to_row:
                    continue
                row = id_to_row[el]
                psi_frame[row, :] += val.data
                ip_count[row]     += 1

            # Average over integration points and store
            for row in range(n_results):
                if ip_count[row] > 0:
                    psi[row, :, m_idx] = psi_frame[row, :] / ip_count[row]

odb.close()
print("ODB closed.")

# ---------------------------------------------------------------------------
# phi: keep full 6-DOF zero-padded layout  [n_nodes*6 x n_modes]
# dof: active DOFs only, single column, format node.dof  (e.g. 1042.2)
# ---------------------------------------------------------------------------
phi_reduced = phi

# Active mask: rows that are non-zero in at least one mode
active_mask = np.any(phi != 0, axis=1)

# Build dof vector: format node_label.dof_index (e.g. node 1042, DOF 2 -> 1042.2)
dof_active = (dof_full[active_mask, 0].astype(np.float32)
              + dof_full[active_mask, 1].astype(np.float32) * 0.1)
dof_active = dof_active.reshape(-1, 1)   # column vector for MATLAB

print("\nphi shape : {} (6 DOFs per node, zeros for inactive DOFs).".format(phi.shape))
print("dof shape : {} (active DOFs only, format = node.dof).".format(dof_active.shape))
print("DOF ordering: sorted by ascending node label, then DOF index 1->6.")
if do_stress:
    print("psi shape : {} (result_locs x 6 components x modes).".format(psi.shape))

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

modes_dict = {
    "phi"   : phi_reduced,
    "fn"    : fn.reshape(-1,1),
    "dof"   : dof_active,
    "nodes" : nodes_mat,
    "elems" : elems_mat,
}

if do_nset:
    modes_dict["nset"]     = nset_mat

if do_stress:
    modes_dict["psi"]      = psi
    modes_dict["elset_id"] = result_ids

savemat(mat_path, modes_dict, do_compression=True)

print("\nVariables in .mat file:")
print("  phi      {} \t- mode shape matrix (6 DOFs/node, zero-padded)".format(phi_reduced.shape))
print("  fn       {} \t\t- natural frequencies (Hz)".format(fn.reshape(-1, 1).shape))
print("  dof      {} \t\t- active DOFs, format node.dof (e.g. 1042.2 = node 1042, U2)".format(dof_active.shape))
print("  nodes    {} \t\t- (label | X | Y | Z)".format(nodes_mat.shape))
print("  elems    {} \t- (label | vtk_type | n1..n8)".format(elems_mat.shape))
if do_stress:
    print("  psi      {} \t- stress tensor (result_locs x 6 components x modes)".format(psi.shape))
    print("  elset_id {} \t\t- (elem_label | integration_point)".format(result_ids.shape))
if do_nset:
    print("  nset     {} \t\t- (label | X | Y | Z)".format(nset_mat.shape))
print("\nElement matrix column layout:")
print("  col 1    : element label")
print("  col 2    : VTK type code  (1=point 3=beam 5=tri 9=quad 10=tet 12=hex ...)")
print("  col 3-10 : node connectivity (zero-padded for elements with < 8 nodes)")
print("\nDOF index legend:  1=U1  2=U2  3=U3  4=UR1  5=UR2  6=UR3")
if do_stress:
    print("psi component order:  1=S11  2=S22  3=S33  4=S12  5=S13  6=S23")
print("\nDone.")
print("\nSaved: {}".format(mat_path))