"""
Find the 3D center of each spine head and assign each center
to one voxel in a regular Cartesian grid.

The grid is shifted to:
    1. Minimize the number of heads sharing a voxel.
    2. Keep head centers close to voxel centers.

The NO point process should use the voxel containing head(0.5).
"""

from pathlib import Path
import csv
import json
import os

import numpy as np
from neuron import h


# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

ROOT = Path(__file__).resolve().parent

VOXEL_SIZE = 0.075

# The grid is tested at SEARCH_STEPS positions along each axis.
# 10 means 10 x 10 x 10 = 1000 possible grid shifts.
SEARCH_STEPS = 10

# Change this if your spine-head section names do not contain "head".
HEAD_NAME_TEXT = "head"


# ------------------------------------------------------------
# Load morphology
# ------------------------------------------------------------

os.chdir(ROOT)

h.load_file("stdrun.hoc")
h.load_file(str(ROOT / "Purk2M0.nrn"))
h.load_file("_generated_PurkMorph101_06.hoc")

head_sections = [
    sec
    for sec in h.allsec()
    if "head" in sec.name().lower()
]

heads_without_3d = [
    sec
    for sec in head_sections
    if int(h.n3d(sec=sec)) == 0
]

print("Heads found:", len(head_sections))
print(
    "Heads with no pt3d before define_shape:",
    len(heads_without_3d),
)


h.define_shape()

sections = list(h.allsec())

print(f"Loaded {len(sections)} total sections")


# ------------------------------------------------------------
# Convert sec(x) into an XYZ coordinate
# ------------------------------------------------------------

def coord_at_x(sec, x):
    """
    Return the interpolated XYZ coordinate at sec(x).

    x is a normalized section location between 0 and 1.
    """

    if not 0.0 <= x <= 1.0:
        raise ValueError(
            f"x must be between 0 and 1, got {x}"
        )

    n = int(h.n3d(sec=sec))

    if n < 2:
        raise ValueError(
            f"{sec.name()} has fewer than 2 3D points"
        )

    total_arc = h.arc3d(n - 1, sec=sec)
    target_arc = x * total_arc

    for i in range(n - 1):
        arc0 = h.arc3d(i, sec=sec)
        arc1 = h.arc3d(i + 1, sec=sec)

        if arc0 <= target_arc <= arc1:

            if arc1 > arc0:
                fraction = (
                    target_arc - arc0
                ) / (
                    arc1 - arc0
                )
            else:
                fraction = 0.0

            x0 = h.x3d(i, sec=sec)
            y0 = h.y3d(i, sec=sec)
            z0 = h.z3d(i, sec=sec)

            x1 = h.x3d(i + 1, sec=sec)
            y1 = h.y3d(i + 1, sec=sec)
            z1 = h.z3d(i + 1, sec=sec)

            return np.array([
                x0 + fraction * (x1 - x0),
                y0 + fraction * (y1 - y0),
                z0 + fraction * (z1 - z0),
            ])

    # Fallback for a tiny floating-point error near x = 1.
    return np.array([
        h.x3d(n - 1, sec=sec),
        h.y3d(n - 1, sec=sec),
        h.z3d(n - 1, sec=sec),
    ])


# ------------------------------------------------------------
# Find the spine-head sections
# ------------------------------------------------------------

def find_spine_heads():
    """
    Find sections whose names contain HEAD_NAME_TEXT.
    """

    heads = []

    for sec in h.allsec():
        if HEAD_NAME_TEXT.lower() in sec.name().lower():
            heads.append(sec)

    return heads


def get_head_centers(head_sections):
    """
    Get the XYZ coordinate of head(0.5) for every spine head.
    """

    centers = []

    for head in head_sections:
        center = coord_at_x(head, 0.5)
        centers.append(center)

    return np.asarray(centers)


# ------------------------------------------------------------
# Voxel-grid functions
# ------------------------------------------------------------

def assign_voxels(coords, origin, voxel_size):
    """
    Return the integer voxel index containing each coordinate.
    """

    return np.floor(
        (coords - origin) / voxel_size
    ).astype(int)


def get_voxel_centers(voxel_indices, origin, voxel_size):
    """
    Return the XYZ center of each assigned voxel.
    """

    return (
        origin
        + (voxel_indices + 0.5) * voxel_size
    )


def count_voxel_collisions(voxel_indices):
    """
    Count how many voxels contain more than one head.
    """

    unique_voxels, counts = np.unique(
        voxel_indices,
        axis=0,
        return_counts=True,
    )

    shared_voxels = int(np.sum(counts > 1))

    extra_heads = int(
        np.sum(counts[counts > 1] - 1)
    )

    return shared_voxels, extra_heads


def find_best_grid(coords, voxel_size, search_steps):
    """
    Try different grid shifts and return the best grid.

    The best grid has:
        1. The fewest heads sharing voxels.
        2. The smallest average distance from heads to voxel centers.
    """

    coordinate_minimum = coords.min(axis=0)

    shifts = np.linspace(
        0.0,
        voxel_size,
        search_steps,
        endpoint=False,
    )

    best_origin = None
    best_voxels = None
    best_extra_heads = None
    best_mean_distance = None

    for shift_x in shifts:
        for shift_y in shifts:
            for shift_z in shifts:

                shift = np.array([
                    shift_x,
                    shift_y,
                    shift_z,
                ])

                origin = coordinate_minimum - shift

                voxels = assign_voxels(
                    coords,
                    origin,
                    voxel_size,
                )

                _, extra_heads = count_voxel_collisions(
                    voxels
                )

                centers = get_voxel_centers(
                    voxels,
                    origin,
                    voxel_size,
                )

                distances = np.linalg.norm(
                    coords - centers,
                    axis=1,
                )

                mean_distance = distances.mean()

                if best_origin is None:
                    better = True

                elif extra_heads < best_extra_heads:
                    better = True

                elif (
                    extra_heads == best_extra_heads
                    and mean_distance < best_mean_distance
                ):
                    better = True

                else:
                    better = False

                if better:
                    best_origin = origin
                    best_voxels = voxels
                    best_extra_heads = extra_heads
                    best_mean_distance = mean_distance

    return best_origin, best_voxels


# ------------------------------------------------------------
# Extract head coordinates
# ------------------------------------------------------------

head_sections = find_spine_heads()

print(f"Found {len(head_sections)} spine-head sections")

if len(head_sections) == 0:
    raise RuntimeError(
        "No spine-head sections were found. "
        "Change HEAD_NAME_TEXT to match your section names."
    )

print("First five head sections:")

for head in head_sections[:5]:
    print(
        head.name(),
        f"L={head.L:.3f}",
        f"diam={head.diam:.3f}",
        f"n3d={int(h.n3d(sec=head))}",
    )

head_coords = get_head_centers(head_sections)

from scipy.spatial import cKDTree

tree = cKDTree(head_coords)

distances, neighbor_indices = tree.query(
    head_coords,
    k=2,
)

nearest_distances = distances[:, 1]

closest_head = np.argmin(nearest_distances)
closest_neighbor = neighbor_indices[closest_head, 1]

print()
print(
    "Minimum distance between head centers:",
    nearest_distances[closest_head],
    "um",
)

print(
    "Closest heads:",
    head_sections[closest_head].name(),
    head_sections[closest_neighbor].name(),
)

print("Coordinates:")
print(head_coords[closest_head])
print(head_coords[closest_neighbor])

rounded_coords = np.round(
    head_coords,
    decimals=6,
)

unique_coords, counts = np.unique(
    rounded_coords,
    axis=0,
    return_counts=True,
)

print("Number of heads:", len(head_coords))
print("Unique head coordinates:", len(unique_coords))
print(
    "Coordinate locations used by multiple heads:",
    np.sum(counts > 1),
)
print(
    "Maximum heads at one coordinate:",
    counts.max(),
)

# ------------------------------------------------------------
# Find and apply the best grid
# ------------------------------------------------------------

grid_origin, head_voxels = find_best_grid(
    head_coords,
    voxel_size=VOXEL_SIZE,
    search_steps=SEARCH_STEPS,
)

head_voxel_centers = get_voxel_centers(
    head_voxels,
    grid_origin,
    VOXEL_SIZE,
)

head_to_voxel_distance = np.linalg.norm(
    head_coords - head_voxel_centers,
    axis=1,
)


# ------------------------------------------------------------
# Check for shared voxels
# ------------------------------------------------------------

unique_voxels, inverse, counts = np.unique(
    head_voxels,
    axis=0,
    return_inverse=True,
    return_counts=True,
)

shared_voxel_indices = np.where(counts > 1)[0]

print()
print(f"Voxel size: {VOXEL_SIZE} um")
print(f"Grid origin: {grid_origin}")
print(f"Number of heads: {len(head_sections)}")
print(f"Number of occupied voxels: {len(unique_voxels)}")
print(
    "Voxels containing multiple heads:",
    len(shared_voxel_indices),
)
print("Maximum heads in one voxel:", counts.max())
print(
    "Mean head-to-voxel-center distance: "
    f"{head_to_voxel_distance.mean():.4f} um"
)
print(
    "Maximum head-to-voxel-center distance: "
    f"{head_to_voxel_distance.max():.4f} um"
)


# Print the heads involved in any voxel collision.
if len(shared_voxel_indices) > 0:

    print("\nShared voxel details:")

    for unique_index in shared_voxel_indices:

        voxel = unique_voxels[unique_index]

        head_indices = np.where(
            inverse == unique_index
        )[0]

        head_names = [
            head_sections[i].name()
            for i in head_indices
        ]

        print(
            f"Voxel {tuple(voxel)} contains "
            f"{len(head_indices)} heads: {head_names}"
        )


# ------------------------------------------------------------
# Save the head-to-voxel mapping
# ------------------------------------------------------------

out_dir = ROOT / "output"
out_dir.mkdir(exist_ok=True)

csv_path = out_dir / "spine_head_voxels.csv"

with csv_path.open("w", newline="") as file:

    writer = csv.writer(file)

    writer.writerow([
        "head_index",
        "head_name",

        "head_x",
        "head_y",
        "head_z",

        "voxel_i",
        "voxel_j",
        "voxel_k",

        "voxel_center_x",
        "voxel_center_y",
        "voxel_center_z",

        "distance_to_voxel_center_um",
    ])

    for i, head in enumerate(head_sections):

        writer.writerow([
            i,
            head.name(),

            head_coords[i, 0],
            head_coords[i, 1],
            head_coords[i, 2],

            head_voxels[i, 0],
            head_voxels[i, 1],
            head_voxels[i, 2],

            head_voxel_centers[i, 0],
            head_voxel_centers[i, 1],
            head_voxel_centers[i, 2],

            head_to_voxel_distance[i],
        ])


# Save the grid settings separately.
grid_path = out_dir / "spine_head_grid.json"

with grid_path.open("w") as file:

    json.dump(
        {
            "voxel_size_um": VOXEL_SIZE,
            "grid_origin_um": grid_origin.tolist(),
            "number_of_heads": len(head_sections),
            "number_of_occupied_voxels": len(unique_voxels),
            "shared_voxels": len(shared_voxel_indices),
            "maximum_heads_in_one_voxel": int(counts.max()),
        },
        file,
        indent=2,
    )


print(f"\nSaved head mapping to: {csv_path}")
print(f"Saved grid information to: {grid_path}")


# ------------------------------------------------------------
# Stop if the grid is not one-head-per-voxel
# ------------------------------------------------------------

if len(shared_voxel_indices) > 0:
    raise RuntimeError(
        "At least two spine-head centers map to the same NO voxel. "
        "Try increasing SEARCH_STEPS or reducing VOXEL_SIZE."
    )

print(
    "\nSuccess: every spine head has one unique NO-producing voxel."
)
