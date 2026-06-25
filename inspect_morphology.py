#!/usr/bin/env python3
"""Load only the HOC morphology from Python and print basic section counts."""
from pathlib import Path
import os
from neuron import h  # , gui
import argparse
import numpy as np

ROOT = Path(__file__).resolve().parent
os.chdir(ROOT)
h.load_file("stdrun.hoc")
h.load_file(str(ROOT / "Purk2M0.nrn"))
h.load_file("_generated_PurkMorph101_06.hoc")
sections = list(h.allsec())
dendrites = [sec for sec in sections if sec.name().startswith("dend")]
print(f"Loaded {len(sections)} sections")
print(f"Dendrites: {len(dendrites)}")
print(f"Soma: L={h.soma.L:g} um, diam={h.soma.diam:g} um")
for sec in sections[:10]:
    print(sec.name(), "L =", sec.L, "diam =", sec.diam, "n3d =", int(h.n3d(sec=sec)))

h.define_shape()


def coord_at_x(sec, x):
    """Interpolate real/defined 3D coordinate along a section at normalized x."""
    n = int(h.n3d(sec=sec))
    if n < 2:
        raise ValueError(f"{sec.name()} has fewer than 2 3D points")

    target_arc = x * h.arc3d(n - 1, sec=sec)

    for i in range(n - 1):
        a0 = h.arc3d(i, sec=sec)
        a1 = h.arc3d(i + 1, sec=sec)

        if a0 <= target_arc <= a1:
            f = (target_arc - a0) / (a1 - a0) if a1 > a0 else 0.0

            return np.array([
                h.x3d(i, sec=sec) + f * (h.x3d(i + 1, sec=sec) - h.x3d(i, sec=sec)),
                h.y3d(i, sec=sec) + f * (h.y3d(i + 1, sec=sec) - h.y3d(i, sec=sec)),
                h.z3d(i, sec=sec) + f * (h.z3d(i + 1, sec=sec) - h.z3d(i, sec=sec)),
            ])

    return np.array([
        h.x3d(n - 1, sec=sec),
        h.y3d(n - 1, sec=sec),
        h.z3d(n - 1, sec=sec),
    ])


def get_spine_anchor_coords():
    """
    Recreate the HOC spine-placement order.

    Template logic:
        forsec sl:
            place spines at x = 0.05, 0.15, ..., 0.95

    Returns:
        coords: (nspines, 3) array
        meta: list of (spine_index, dend_name, dend_x)
    """
    spine_xs = np.arange(0.05, 1.0, 0.1)

    coords = []
    meta = []
    spine_i = 0

    for dend in h.sl:
        for x in spine_xs:
            coord = coord_at_x(dend, float(x))
            coords.append(coord)
            meta.append((spine_i, dend.name(), float(x)))
            spine_i += 1

    return np.asarray(coords), meta


def compare_to_regular_grid(coords, grid_spacing):
    """
    Find nearest point on a regular Cartesian grid and distance to it.

    Grid origin is coords.min(axis=0), so offsets are relative to the
    morphology bounding box.
    """
    mins = coords.min(axis=0)

    grid_indices = np.round((coords - mins) / grid_spacing)
    nearest_grid = mins + grid_indices * grid_spacing

    offsets = coords - nearest_grid
    distances = np.linalg.norm(offsets, axis=1)

    return nearest_grid, offsets, distances

parser = argparse.ArgumentParser()
parser.add_argument("--grid-spacing", type=float, default=2.0)
args = parser.parse_args()

spine_coords, spine_meta = get_spine_anchor_coords()

print(f"Spine anchor coords: {spine_coords.shape}")
print("First 10 spine anchors:")
for i in range(min(10, len(spine_meta))):
    spine_i, dend_name, dend_x = spine_meta[i]
    x, y, z = spine_coords[i]
    print(f"{spine_i}: {dend_name}({dend_x:.2f}) -> ({x:.3f}, {y:.3f}, {z:.3f})")

nearest_grid, offsets, distances = compare_to_regular_grid(
    spine_coords,
    args.grid_spacing,
)

print(f"\nGrid spacing: {args.grid_spacing} um")
print(f"Mean offset: {distances.mean():.4f} um")
print(f"Median offset: {np.median(distances):.4f} um")
print(f"Max offset: {distances.max():.4f} um")
print(f"95th percentile offset: {np.percentile(distances, 95):.4f} um")

out_dir = ROOT / "output"
out_dir.mkdir(exist_ok=True)

np.savetxt(
    out_dir / "spine_anchor_coords.csv",
    spine_coords,
    delimiter=",",
    header="x,y,z",
    comments="",
)

np.savetxt(
    out_dir / f"spine_grid_offsets_{args.grid_spacing:g}um.csv",
    np.column_stack([spine_coords, nearest_grid, offsets, distances]),
    delimiter=",",
    header="spine_x,spine_y,spine_z,grid_x,grid_y,grid_z,dx,dy,dz,distance_um",
    comments="",
)

#
# no_voxels = []
# coords = []
# for sec in h.allsec():
#     if 'spine_head' in sec.name():
#         vox = h.no_voxel(sec(0.5))
#         no_voxels.append(vox)


# shape = h.Shape()

