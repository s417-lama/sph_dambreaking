#!/usr/bin/env python3

import sys

scale = int(sys.argv[1]) if len(sys.argv) > 1 else 1

l0      = 0.55 / 30.0 / scale
n_wall  = 2

fluid_x = 67 * scale
fluid_y = 54 * scale
fluid_z = 30 * scale

world_x = 175 * scale
world_y = 54 * scale
world_z = 55 * scale

box_x = 10 * scale
box_y = 20 * scale
box_z = 9 * scale

box_offset_x = 130 * scale

# fluid
for x in range(fluid_x):
    for y in range(fluid_y):
        for z in range(fluid_z):
            print((x + 1) * l0, (y + 1) * l0, (z + 1) * l0, 1)

# bottom wall
for x in range(-n_wall, world_x + n_wall):
    for y in range(-n_wall, world_y + n_wall):
        for z in range(-n_wall, 0):
            print((x + 1) * l0, (y + 1) * l0, (z + 1) * l0, 2)

# xz walls
for x in range(-n_wall, world_x + n_wall):
    for y in range(-n_wall, 0):
        for z in range(world_z):
            print((x + 1) * l0, (y + 1) * l0, (z + 1) * l0, 2)

for x in range(-n_wall, world_x + n_wall):
    for y in range(world_y, world_y + n_wall):
        for z in range(world_z):
            print((x + 1) * l0, (y + 1) * l0, (z + 1) * l0, 2)

# yz walls
for x in range(-n_wall, 0):
    for y in range(world_y):
        for z in range(world_z):
            print((x + 1) * l0, (y + 1) * l0, (z + 1) * l0, 2)

for x in range(world_x, world_x + n_wall):
    for y in range(world_y):
        for z in range(world_z):
            print((x + 1) * l0, (y + 1) * l0, (z + 1) * l0, 2)

# box
for x in range(box_offset_x, box_offset_x + box_x):
    for y in range(world_y // 2 - box_y // 2, world_y // 2 + box_y // 2):
        for z in range(box_z):
            print((x + 1) * l0, (y + 1) * l0, (z + 1) * l0, 2)
