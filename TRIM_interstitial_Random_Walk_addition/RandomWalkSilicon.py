import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D
import random
from itertools import product
from math import floor

random.seed(0)
np.random.seed(0)

def plot_lattice(displacements, oxygen_map, divacancies, vac_oxygens, trivavancies, step):
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    if displacements:
        dx, dy, dz = zip(*displacements)
        ax.scatter(dx, dy, dz, c='blue', marker='o', label=str(len(displacements))+' Displacements')
    if oxygen_map:
        dx, dy, dz = zip(*oxygen_map)
        ax.scatter(dx, dy, dz, c='yellow', marker='.', label=f"{len(oxygen_map)} Oxygens")
    if divacancies:
        dvx, dvy, dvz = zip(*divacancies)
        ax.scatter(dvx, dvy, dvz, c='red', marker='x', label=str(len(divacancies))+' Divacancies')
    if vac_oxygens:
        vox, voy, voz = zip(*vac_oxygens)
        ax.scatter(vox, voy, voz, c='green', marker='^', label=str(len(vac_oxygens))+' VO')
    if trivavancies:
        dvx, dvy, dvz = zip(*trivavancies)
        ax.scatter(dvx, dvy, dvz, c='purple', marker='+', label=f"{len(trivavancies)} Trivacancies")
    # Set the limits for x, y, z axes
    ax.set_xlim(0, lattice_size - 1)
    ax.set_ylim(0, lattice_size - 1)
    ax.set_zlim(0, lattice_size - 1)
    
    # Set tick labels to show only whole numbers
    ticks = range(0, lattice_size, 1)  # Change the step if necessary for very large lattices
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_zticks(ticks)

    # Axis labels and legend
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.legend()

    # Title and saving the plot
    plt.title(f'Step {step}')
    plt.savefig(f'{output_dir}/step_{step}.png')
    plt.close()

def calculate_neighbors(position, max_distance):
    x, y, z = position
    neighbors = []
    for dx, dy, dz in product(range(-max_distance, max_distance + 1), repeat=3):  # Range based on max_distance
        if dx == 0 and dy == 0 and dz == 0:
            continue  # Skip the current position itself
        if dx**2 + dy**2 + dz**2 <= max_distance**2:
            new_position = (x + dx, y + dy, z + dz)
            if 0 <= new_position[0] < lattice_size and 0 <= new_position[1] < lattice_size and 0 <= new_position[2] < lattice_size:
                neighbors.append(new_position)
    return neighbors
# Create directory for output figures if it does not exist
output_dir = './DemoFig2'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

################PROP SI DEFINITION###################
# Define the Silicon lattice
# a = 5.431  # lattice parameter of Silicon in Ångstroms
# lattice = Lattice.cubic(a)
# structure = Structure(lattice, ["Si", "Si"], [[0, 0, 0], [0.25, 0.25, 0.25]], coords_are_cartesian=False)
# # Assuming a certain size of the simulation cell (number of unit cells in each dimension)
# cell_size = 10  # 10x10x10 unit cells
# super_lattice = structure * (cell_size, cell_size, cell_size)  # Create a larger lattice

# # For demonstration, you can consider all Si positions initially as potential displacement positions
# displacements = {tuple(map(int, site.frac_coords * a * cell_size)) for site in super_lattice}

# def move_displacement(displacement):
#     # Example simple move, refine based on actual physical movement rules
#     moves = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, 0, 0), (0, -1, 0), (0, 0, -1)]
#     move = random.choice(moves)
#     return tuple(displacement[i] + move[i] for i in range(3))
# ############CLASSIC 3D DEFINITION###################
# Lattice size and initialization
lattice_size = 100  # 20x20x20 lattice
GRID_UNIT = 0.235 # nm
STEPS = 200  # Number of steps in the simulation

PROB_DIVACANCY = 0.107
DIVACANCY_NEIGHBOURHOOD = 3
PROB_VAC_OXYGEN = 0.029
VAC_OXYGEN_NEIGHBOURHOOD = 2
PROB_TRIVACANCY = 0.226
TRIVACANCY_NEIGHBOURHOOD = 4

#Initial number of displacements
num_displacements = 500
displacements = set()
while len(displacements) < num_displacements:
    displacements.add(tuple(np.random.randint(0, lattice_size, size=3)))
#Initial number of oxygen
SAMPLE_SIZE_X = 26000000  # nm
SAMPLE_SIZE_Y = 26000000  # nm
SAMPLE_SIZE_Z = 150000  # nm
SAMPLE_OXYGEN_COUNT = 10**18

oxygen_map = set()
number_of_oxygens = floor(SAMPLE_OXYGEN_COUNT / (SAMPLE_SIZE_X * SAMPLE_SIZE_Y * SAMPLE_SIZE_Z) * (lattice_size * GRID_UNIT)**3)
while len(oxygen_map) < number_of_oxygens:
    oxygen_map.add(tuple(np.random.randint(0, lattice_size, size=3)))

divacancies = set()  # Initially no divacancies
vac_oxygens = set()
trivacancies = set()
displacement_count = []
divacancy_count = []
vac_oxygens_count = []
trivacancy_count = []

for step in range(STEPS):
    plot_lattice(displacements, oxygen_map, divacancies, vac_oxygens, trivacancies, step)

    next_divacancies = set(divacancies)
    next_vac_oxygens = set(vac_oxygens)
    next_trivacancies = set(trivacancies)

    # Check for potential divacancies
    remaining_displacements = set()
    while displacements:
        displacement = displacements.pop()

        neighbours_for_divacancy = calculate_neighbors(displacement, DIVACANCY_NEIGHBOURHOOD)  
        neighbours_for_oxygens = calculate_neighbors(displacement, VAC_OXYGEN_NEIGHBOURHOOD)
        neighbours_for_trivacancy = calculate_neighbors(displacement, TRIVACANCY_NEIGHBOURHOOD)

        interaction_flag = False

        #neighbors = [tuple(np.array(d) + np.array(delta)) for delta in [(1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)]]
        if (neighbouring_divacancies := [n for n in neighbours_for_trivacancy if n in divacancies]):
            for divacancy in neighbouring_divacancies:
                if random.random() < PROB_TRIVACANCY:
                    next_trivacancies.add(divacancy)
                    next_divacancies.discard(divacancy)
                    interaction_flag = True
                    break

        elif (neighbouring_displacements := [n for n in neighbours_for_divacancy if n in displacements]):
            for neighbouring_displacement in neighbouring_displacements:
                if random.random() < PROB_DIVACANCY:
                    next_divacancies.add(displacement)  # this could be added between displacement and neighbouring_displacement
                    displacements.discard(neighbouring_displacement)
                    interaction_flag = True
                    break

        # is priority of divacancy higher than vac oxygen?
        elif (neighbouring_oxygens := [n for n in neighbours_for_oxygens if n in oxygen_map]):
            for oxygen in neighbouring_oxygens:
                if random.random() < PROB_VAC_OXYGEN:
                    next_vac_oxygens.add(oxygen)
                    oxygen_map.discard(oxygen)
                    interaction_flag = True
                    break

        if not interaction_flag:
            remaining_displacements.add(displacement)

    # update sets
    divacancies = next_divacancies
    vac_oxygens = next_vac_oxygens
    trivacancies = next_trivacancies

    # Perform random walks
    moved_displacements = set()
    for displacement in remaining_displacements:
        move = np.random.randint(-1, 2, size=3)
        new_position = tuple(np.array(displacement) + move)
        if 0 <= min(new_position) < lattice_size and max(new_position) < lattice_size:
            if new_position not in divacancies and new_position not in vac_oxygens and new_position not in moved_displacements:
                moved_displacements.add(new_position)
            else:
                moved_displacements.add(displacement)
        # Boundary or invalid move cases are ignored (displacement lost)
    displacements = moved_displacements

    # Record counts
    displacement_count.append(len(displacements))
    divacancy_count.append(len(divacancies))
    vac_oxygens_count.append(len(vac_oxygens))
    trivacancy_count.append(len(trivacancies))

plot_lattice(displacements, oxygen_map, divacancies, vac_oxygens, trivacancies, STEPS)

# Final plot of divacancy and displacement trends
plt.plot(range(STEPS), displacement_count, label='Displacements')
plt.plot(range(STEPS), divacancy_count, label='Divacancies')
plt.plot(range(STEPS), vac_oxygens_count, label='VO')
plt.plot(range(STEPS), trivacancy_count, label='Trivacancies')
plt.xlabel('Step')
plt.ylabel('Count')
plt.title('Simulation Results Over Time')
plt.legend()
plt.savefig(f'{output_dir}/final_counts.png')
plt.show()
