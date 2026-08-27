#!/usr/bin/env python3
"""Small, standalone NO-voxel diffusion experiment.

This file does not change the Purkinje model.  It loads the existing model,
builds a local Cartesian grid around one spine head, maps nearby heads into
that grid, and records the NO concentration in their voxels.

The numerical values for D, decay, and production are starting values only;
replace them with the values appropriate for the experiment.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
from neuron import h

import run_purk_pr


ROOT = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# Experiment settings -- edit these values, then run this file in your editor.
# ---------------------------------------------------------------------------

PATTERN_INDEX = 0
SOURCE_SPINE = 159

VOXEL_SIZE = 0.3          # um
GRID_HALF_WIDTH = 3.0     # um; local cube extends this far from SOURCE_SPINE

DIFFUSION = 3.3           # um^2/ms; replace with the desired NO value
DECAY_RATE = 0.01         # 1/ms
PRODUCTION_RATE = 1.0     # nM/ms per pattern-active spine

SOURCE_START = 10.0       # ms
SOURCE_STOP = 20.0        # ms
TSTOP = 50.0              # ms
DT = 0.005                # ms


def coordinate_at_x(section, position=0.5):
    """Return the interpolated 3-D coordinate of section(position)."""
    count = int(h.n3d(sec=section))
    if count < 2:
        raise ValueError(f"{section.name()} has fewer than two 3-D points")

    total_length = float(h.arc3d(count - 1, sec=section))
    target = position * total_length

    for point in range(count - 1):
        arc0 = float(h.arc3d(point, sec=section))
        arc1 = float(h.arc3d(point + 1, sec=section))
        if arc0 <= target <= arc1:
            fraction = 0.0 if arc1 == arc0 else (target - arc0) / (arc1 - arc0)
            xyz0 = np.array([
                h.x3d(point, sec=section),
                h.y3d(point, sec=section),
                h.z3d(point, sec=section),
            ])
            xyz1 = np.array([
                h.x3d(point + 1, sec=section),
                h.y3d(point + 1, sec=section),
                h.z3d(point + 1, sec=section),
            ])
            return xyz0 + fraction * (xyz1 - xyz0)

    return np.array([
        h.x3d(count - 1, sec=section),
        h.y3d(count - 1, sec=section),
        h.z3d(count - 1, sec=section),
    ])


def get_spine_heads_and_coordinates():
    """Return spine-head sections in numeric array order and their centers."""
    h.define_shape()
    heads = [
        section for section in h.allsec()
        if section.name().startswith("spine_head[")
    ]
    heads.sort(key=lambda section: int(section.name().split("[")[1].split("]")[0]))
    coordinates = np.asarray([coordinate_at_x(section) for section in heads])
    return heads, coordinates


def build_local_grid(center, half_width, voxel_size, diffusion, decay):
    """Build and wire a cubic grid centered around one spine head."""
    origin = np.floor((center - half_width) / voxel_size) * voxel_size
    upper = np.ceil((center + half_width) / voxel_size) * voxel_size
    shape = np.ceil((upper - origin) / voxel_size).astype(int)
    total_voxels = int(np.prod(shape))

    print(f"Grid origin: {origin}")
    print(f"Grid shape: {tuple(shape)}")
    print(f"Grid voxels: {total_voxels}")

    # POINT_PROCESS objects must be hosted by a NEURON section.  They all may
    # share this artificial section because their spatial relationships are
    # defined by the POINTER connections below, not by cable geometry.
    host = h.Section(name="no_grid_host")
    host.L = 1.0
    host.diam = 1.0

    voxels = {}
    for i in range(shape[0]):
        for j in range(shape[1]):
            for k in range(shape[2]):
                voxel = h.no_voxel(host(0.5))
                voxel.conc0 = 0.0
                voxel.lam = decay
                voxel.F = 0.0
                voxels[(i, j, k)] = voxel

    neighbor_fields = [
        ((1, 0, 0), "conc_xp", "dx_pos"),
        ((-1, 0, 0), "conc_xn", "dx_neg"),
        ((0, 1, 0), "conc_yp", "dy_pos"),
        ((0, -1, 0), "conc_yn", "dy_neg"),
        ((0, 0, 1), "conc_zp", "dz_pos"),
        ((0, 0, -1), "conc_zn", "dz_neg"),
    ]
    diffusion_rate = diffusion / voxel_size**2

    for (i, j, k), voxel in voxels.items():
        for (di, dj, dk), pointer_name, rate_name in neighbor_fields:
            neighbor = voxels.get((i + di, j + dj, k + dk))
            if neighbor is None:
                # Reflecting (no-flux) outer boundary.
                neighbor = voxel
                rate = 0.0
            else:
                rate = diffusion_rate

            h.setpointer(neighbor._ref_conc, pointer_name, voxel)
            setattr(voxel, rate_name, rate)

    return host, voxels, origin, shape


def voxel_index(coordinate, origin, voxel_size):
    return tuple(np.floor((coordinate - origin) / voxel_size).astype(int))


def map_heads_in_grid(coordinates, origin, shape, voxel_size):
    """Map heads inside the local grid; heads outside it are omitted."""
    mapping = {}
    for spine_index, coordinate in enumerate(coordinates):
        index = voxel_index(coordinate, origin, voxel_size)
        if all(0 <= index[axis] < shape[axis] for axis in range(3)):
            mapping[spine_index] = index
    return mapping


def configure_sources(voxels, mapping, rate_per_spine, start, stop):
    """Sum pattern-selected spine sources and play them into voxel.F."""
    sources_by_voxel = defaultdict(float)
    for spine_index, index in mapping.items():
        if h.p.x[spine_index] != 0:
            sources_by_voxel[index] += rate_per_spine

    # Vector.play objects must remain alive throughout the simulation.
    players = []
    for index, total_rate in sources_by_voxel.items():
        times = h.Vector([0.0, start, stop])
        values = h.Vector([0.0, total_rate, 0.0])
        values.play(voxels[index]._ref_F, times)
        players.append((times, values))

    return players, sources_by_voxel


def main():
    run_purk_pr.load_model(PATTERN_INDEX)
    heads, coordinates = get_spine_heads_and_coordinates()

    if not 0 <= SOURCE_SPINE < len(heads):
        raise ValueError(f"SOURCE_SPINE must be between 0 and {len(heads) - 1}")

    source_coordinate = coordinates[SOURCE_SPINE]
    host, voxels, origin, shape = build_local_grid(
        source_coordinate,
        GRID_HALF_WIDTH,
        VOXEL_SIZE,
        DIFFUSION,
        DECAY_RATE,
    )
    mapping = map_heads_in_grid(coordinates, origin, shape, VOXEL_SIZE)

    heads_by_voxel = defaultdict(list)
    for spine_index, index in mapping.items():
        heads_by_voxel[index].append(spine_index)

    shared = {index: members for index, members in heads_by_voxel.items() if len(members) > 1}
    print(f"Spines inside local grid: {len(mapping)}")
    print(f"Occupied local voxels: {len(heads_by_voxel)}")
    print(f"Shared local voxels: {len(shared)}")

    players, source_totals = configure_sources(
        voxels,
        mapping,
        PRODUCTION_RATE,
        SOURCE_START,
        SOURCE_STOP,
    )
    print(f"Pattern-active source voxels: {len(source_totals)}")

    time = h.Vector().record(h._ref_t)
    recordings = {}
    for spine_index, index in mapping.items():
        recordings[spine_index] = h.Vector().record(voxels[index]._ref_conc)

    h.dt = DT
    h.tstop = TSTOP
    h.run()

    output = ROOT / "output"
    output.mkdir(exist_ok=True)
    path = output / f"voxel_test_spine{SOURCE_SPINE}_pattern{PATTERN_INDEX}.csv"

    ordered_spines = sorted(mapping)
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["t_ms"] + [f"spine_{index}_NO_nM" for index in ordered_spines])
        for sample in range(int(time.size())):
            writer.writerow(
                [time.x[sample]]
                + [recordings[index].x[sample] for index in ordered_spines]
            )

    # Keep these objects referenced until after h.run().
    _keep_alive = (host, voxels, players, recordings)
    print(f"Saved: {path}")


if __name__ == "__main__":
    main()
