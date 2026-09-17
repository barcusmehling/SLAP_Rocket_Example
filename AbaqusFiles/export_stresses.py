"""
odb_stress_to_csv.py
--------------------
Extracts modal stress results from an Abaqus ODB file for a specified
element set and saves them to a CSV file matching the Abaqus Report tab
format:
  ODB Name, Step, Frame, Part Instance Name, Element Label, IntPt,
  X, Y, Z, Section Name, Material Name, Section Point,
  S-S11, S-S22, S-S33, S-S12, S-S13, S-S23

Usage (must be run inside Abaqus Python):
  abaqus python odb_stress_to_csv.py --odb path/to/result.odb --elset NAME [options]

Options:
  --step     Step name            (default: last step)
  --instance Assembly instance    (default: first instance)
  --csv      Output .csv path     (default: <odb_name>_Stress.csv)

Dependencies:
  odbAccess  - bundled with Abaqus
  numpy      - bundled with Abaqus Python
"""

import sys
import os
import csv
import argparse
import numpy as np

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Export Abaqus modal stress ODB to CSV")
parser.add_argument("--odb",      required=True, help="Path to the .odb file")
parser.add_argument("--step",     default=None,  help="Step name (default: last step)")
parser.add_argument("--elset",    required=True, help="Element set name to extract stress from")
parser.add_argument("--mat",      default=None,  help="Output .mat file")
parser.add_argument("--instance", default=None,  help="Assembly instance name")
args = parser.parse_args()

odb_path    = os.path.abspath(args.odb)
odb_name    = os.path.splitext(os.path.basename(odb_path))[0]
default_dir = os.path.abspath(os.path.join(os.path.dirname(odb_path), "..", "ModeShapes"))
if not os.path.exists(default_dir):
    os.makedirs(default_dir)
mat_path = args.mat or os.path.join(default_dir, odb_name + "_Stress.mat")
elset_name  = args.elset

# ---------------------------------------------------------------------------
# Open ODB
# ---------------------------------------------------------------------------
try:
    from odbAccess import openOdb
except ImportError:
    sys.exit(
        "ERROR: odbAccess not found.\n"
        "Run this script with:  abaqus python odb_stress_to_csv.py --odb <file.odb> --elset <name>"
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
# Select instance(s) and resolve element set
# ---------------------------------------------------------------------------
instance_names = list(odb.rootAssembly.instances.keys())

# Try to find the elset: first on the specified/default instance,
# then at the assembly level (spans multiple instances).
elset_region  = None
instance      = None

if args.instance:
    # User specified an instance — use it directly
    if args.instance not in odb.rootAssembly.instances:
        odb.close()
        sys.exit("ERROR: Instance '{}' not found. Available: {}".format(
            args.instance, instance_names))
    instance     = odb.rootAssembly.instances[args.instance]
    if elset_name not in instance.elementSets:
        odb.close()
        sys.exit("ERROR: Element set '{}' not found on instance '{}'.\n"
                 "       If it is an assembly-level set, omit --instance.".format(
                 elset_name, args.instance))
    elset_region = instance.elementSets[elset_name]
    print("Using instance: '{}'".format(instance.name))

else:
    # No instance specified — check assembly-level sets first
    if elset_name in odb.rootAssembly.elementSets:
        elset_region = odb.rootAssembly.elementSets[elset_name]
        print("Element set '{}' found at assembly level (may span multiple instances).".format(
            elset_name))
        # instance stays None — we will resolve per-element below

    else:
        # Fall back: search each instance individually
        for iname in instance_names:
            inst = odb.rootAssembly.instances[iname]
            if elset_name in inst.elementSets:
                instance     = inst
                elset_region = inst.elementSets[elset_name]
                print("Element set '{}' found on instance '{}'.".format(elset_name, iname))
                break

        if elset_region is None:
            odb.close()
            sys.exit("ERROR: Element set '{}' not found on any instance or at assembly level.\n"
                     "       Available instances: {}".format(elset_name, instance_names))

# ---------------------------------------------------------------------------
# Build element label -> instance lookup (needed for assembly-level elsets)
# ---------------------------------------------------------------------------
# If instance is known, this is trivial. If the elset spans multiple instances,
# we map each element label to its owning instance for coordinate/property lookups.
elem_to_instance = {}
elset_labels     = set()

for el in elset_region.elements:
    elset_labels.add(el.label)
    # Assembly-level elements carry an instanceName attribute
    if instance is None:
        iname = el.instanceName
        elem_to_instance[el.label] = odb.rootAssembly.instances[iname]
    else:
        elem_to_instance[el.label] = instance

print("Element set '{}': {} elements across {} instance(s).".format(
    elset_name, len(elset_labels),
    len(set(i.name for i in elem_to_instance.values()))))

# Deduplicate instances by name
seen  = set()
insts = []
for inst in elem_to_instance.values():
    if inst.name not in seen:
        seen.add(inst.name)
        insts.append(inst)

# ---------------------------------------------------------------------------
# Validate element set and build label lookup
# ---------------------------------------------------------------------------
elset_name = args.elset
if elset_name not in instance.elementSets:
    odb.close()
    sys.exit("ERROR: Element set '{}' not found. Available: {}".format(
        elset_name, list(instance.elementSets.keys())))

elset_labels = set(el.label for el in instance.elementSets[elset_name].elements)
print("Element set '{}': {} elements.".format(elset_name, len(elset_labels)))

# ---------------------------------------------------------------------------
# Build node coordinate lookup  — per instance
# ---------------------------------------------------------------------------
print("Building node coordinate lookup ...")

# Keyed by (instance_name, node_label) to avoid collisions across instances
node_coord = {}
for inst in insts:
    for node in inst.nodes:
        coords = list(node.coordinates)
        while len(coords) < 3:
            coords.append(0.0)
        node_coord[(inst.name, node.label)] = coords

# ---------------------------------------------------------------------------
# Build element -> nodes / section / material lookup
# ---------------------------------------------------------------------------
print("Building element property lookup ...")
elem_nodes    = {}
elem_section  = {}
elem_material = {}

for inst in insts:
    for elem in inst.elements:
        if elem.label in elset_labels:
            elem_nodes[elem.label] = (inst.name, list(elem.connectivity))

    for sa in inst.sectionAssignments:
        try:
            sec   = sa.section
            sname = sec.name
            try:
                mname = sec.material
            except AttributeError:
                mname = ""
        except AttributeError:
            # Some Abaqus versions expose sectionName directly on the assignment
            try:
                sname = sa.sectionName
                mname = ""
            except AttributeError:
                sname = ""
                mname = ""

        try:
            region_elems = sa.region.elements
        except AttributeError:
            continue

        for el in region_elems:
            if el.label in elset_labels:
                elem_section[el.label]  = sname
                elem_material[el.label] = mname

# ---------------------------------------------------------------------------
# First-pass: discover all (element_label, integration_point) pairs
# ---------------------------------------------------------------------------
print("Discovering result locations from first mode frame ...")

frames      = step.frames
mode_frames = [f for f in frames if f.frameValue > 0]
n_modes     = len(mode_frames)
print("  {} modes found.".format(n_modes))

first_frame  = mode_frames[0]
elset_region = instance.elementSets[elset_name]
subset       = first_frame.fieldOutputs["S"].getSubset(region=elset_region)

result_ids = []
for val in subset.values:
    result_ids.append((val.elementLabel, val.integrationPoint))

result_ids.sort(key=lambda x: (x[0], x[1]))
n_results = len(result_ids)
print("  {} result locations (element/IP pairs) found.".format(n_results))

# ---------------------------------------------------------------------------
# Extract stresses
# psi : [n_results x 6 x n_modes]  (element/IP x stress component x mode)
#   component order: S11 S22 S33 S12 S13 S23
# ---------------------------------------------------------------------------
print("Extracting stresses ...")

n_results = len(result_ids)
psi       = np.zeros((n_results, 6, n_modes), dtype=np.float64)

for m_idx, frame in enumerate(mode_frames):

    # --- Parse frequency from frame description ----------------------------
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
    print("  Mode {:>3}: {:>12.4f} Hz  | {}".format(m_idx + 1, freq_hz, desc.strip()))

    # --- Extract stress for this frame -------------------------------------
    if "S" not in frame.fieldOutputs:
        print("    WARNING: 'S' field not found in mode {}, skipping.".format(m_idx + 1))
        continue

    subset = frame.fieldOutputs["S"].getSubset(region=elset_region)

    stress_map = {}
    for val in subset.values:
        key = (int(val.elementLabel), int(val.integrationPoint))
        stress_map[key] = val.data    # (S11, S22, S33, S12, S13, S23)

    for r_idx, (el_label, ip) in enumerate(result_ids):
        data = stress_map.get((el_label, ip))
        if data is None:
            continue
        psi[r_idx, :, m_idx] = data

odb.close()
print("ODB closed.")

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

mat_path = args.csv.replace(".csv", ".mat") if args.csv else os.path.join(
    default_dir, odb_name + "_Stress.mat")

save_dict = {
    "psi"     : psi,                    # [n_results x 6 x n_modes]
    "fn"      : fn.reshape(-1, 1),      # [n_modes x 1]  column vector in MATLAB
    "elset_id": result_ids,             # [n_results x 2]  (elem_label | int_point)
}

savemat(mat_path, save_dict, do_compression=True)

print("\nVariables in .mat file:")
print("  psi      {} - stress tensor (result_locs x 6 components x modes)".format(psi.shape))
print("  fn       {} - natural frequencies (Hz)".format(fn.reshape(-1, 1).shape))
print("  elset_id {} - (elem_label | integration_point)".format(result_ids.shape))
print("\nComponent order:  1=S11  2=S22  3=S33  4=S12  5=S13  6=S23")
print("\nDone.")
print("\nSaved: {}".format(mat_path))