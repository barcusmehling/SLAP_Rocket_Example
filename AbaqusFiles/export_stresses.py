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
parser.add_argument("--csv",      default=None,  help="Output .csv file")
parser.add_argument("--instance", default=None,  help="Assembly instance name")
args = parser.parse_args()

odb_path    = os.path.abspath(args.odb)
odb_name    = os.path.splitext(os.path.basename(odb_path))[0]
default_dir = os.path.abspath(os.path.join(os.path.dirname(odb_path), "..", "ModeShapes"))
if not os.path.exists(default_dir):
    os.makedirs(default_dir)
csv_path = args.csv or os.path.join(default_dir, odb_name + "_Stress.csv")

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
# Build node coordinate lookup  node_label -> (X, Y, Z)
# ---------------------------------------------------------------------------
print("Building node coordinate lookup ...")
node_coord = {}
for node in instance.nodes:
    coords = list(node.coordinates)
    while len(coords) < 3:
        coords.append(0.0)
    node_coord[node.label] = coords

# ---------------------------------------------------------------------------
# Build element -> (section_name, material_name, node_labels) lookup
# ---------------------------------------------------------------------------
print("Building element property lookup ...")
elem_nodes = {}
for elem in instance.elements:
    if elem.label in elset_labels:
        elem_nodes[elem.label] = list(elem.connectivity)

# Section/material: walk sectionAssignments on the instance
elem_section  = {}   # elem_label -> section name
elem_material = {}   # elem_label -> material name
for sa in instance.sectionAssignments:
    sec  = sa.section
    sname = sec.name
    try:
        mname = sec.material
    except AttributeError:
        mname = ""
    region = sa.region
    for el in region.elements:
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
# Extract stresses and write CSV
# ---------------------------------------------------------------------------
print("Extracting stresses and writing CSV ...")

CSV_HEADER = [
    "ODB Name", "Step", "Frame", "Part Instance Name",
    " Element Label", "         IntPt",
    "X", "Y", "Z",
    "Section Name", "Material Name", "Section Point",
    "         S-S11", "         S-S22", "         S-S33",
    "         S-S12", "         S-S13", "         S-S23",
]

with open(csv_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(CSV_HEADER)

    for m_idx, frame in enumerate(mode_frames):

        # --- Parse frequency from frame description ------------------------
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

        frame_str = desc.strip()   # e.g. "Mode  1: Value = ... Freq = ... (cycles/time)"
        print("  Mode {:>3}: {:>12.4f} Hz  | {}".format(m_idx + 1, freq_hz, frame_str))

        # --- Extract stress for this frame ---------------------------------
        if "S" not in frame.fieldOutputs:
            print("    WARNING: 'S' field not found in mode {}, skipping.".format(m_idx + 1))
            continue

        subset = frame.fieldOutputs["S"].getSubset(region=elset_region)

        # Collect values into a dict keyed by (elem, ip) for ordered writing
        stress_map = {}
        for val in subset.values:
            key = (int(val.elementLabel), int(val.integrationPoint))
            stress_map[key] = val

        # Write rows in sorted (elem_label, int_point) order
        for (el_label, ip) in result_ids:
            key = (el_label, ip)
            if key not in stress_map:
                continue

            val   = stress_map[key]
            data  = val.data          # (S11, S22, S33, S12, S13, S23)

            # Node coordinates: average over element connectivity for the IP centroid
            # (Abaqus does not expose IP coordinates directly; use element node average)
            conn   = elem_nodes.get(el_label, [])
            coords = [node_coord.get(n, [0.0, 0.0, 0.0]) for n in conn]
            if coords:
                x = sum(c[0] for c in coords) / len(coords)
                y = sum(c[1] for c in coords) / len(coords)
                z = sum(c[2] for c in coords) / len(coords)
            else:
                x = y = z = 0.0

            sname = elem_section.get(el_label,  "")
            mname = elem_material.get(el_label, "")

            writer.writerow([
                odb_path,
                step.name,
                frame_str,
                '"{}"'.format(instance.name),
                el_label,
                ip,
                x, y, z,
                '"{}"'.format(sname),
                '"{}"'.format(mname),
                '""',          # Section Point — empty as in the reference
                data[0], data[1], data[2],
                data[3], data[4], data[5],
            ])

odb.close()
print("ODB closed.")

print("\nDone.")
print("\nSaved: ../ModeShapes/{}_Stress.csv".format(odb_name))