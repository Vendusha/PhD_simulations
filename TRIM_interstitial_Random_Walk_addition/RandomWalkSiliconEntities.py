from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D
import random
from itertools import product
from math import floor, ceil
from enum import Enum
from dataclasses import dataclass
from random import shuffle, random, seed, choice
from tqdm import tqdm

seed(0)
np.random.seed(0)


# Lattice size and initialization
LATTICE_SIZE = 100  # 23.5 x 23.5x 23.5 lattice
GRID_UNIT = 0.235 # nm
STEPS = 200  # Number of steps in the simulation


class EntityKind(Enum):
    VACANCY = 1
    INTERSTITIAL = 2
    DIVACANCY = 3
    TRIVACANCY = 4
    OXYGEN = 5
    VAC_OXYGEN = 6
    DIVAC_OXYGEN = 7
    TRIVAC_OXYGEN = 8
    DIINTERSTITIAL = 9
    CARBON = 10
    CARBON_INTERSTITIAL = 11

MOVING_ENTITIES = {EntityKind.VACANCY, EntityKind.INTERSTITIAL}

class EntitiesMeta():
    container = None

    def __init__(self):
        for entity_kind in EntityKind:
            setattr(self, entity_kind.name.lower(), self.container())
    
    def get_kind(self, entity_kind: EntityKind):
        return getattr(self, entity_kind.name.lower())


class Entities(EntitiesMeta):
    container = set

    def copy(self):
        entities = Entities()
        for attr, value in vars(self).items():
            setattr(entities, attr, set(value))
        return entities

    def report_counts(self):
        return {attr: len(value) for attr, value in vars(self).items()}


@dataclass
class InteractionKind():
    interactant_a_kind: EntityKind
    interactant_b_kind: EntityKind
    radius: int
    probability: float
    process: callable

    def __post_init__(self):
        if self.interactant_a_kind not in MOVING_ENTITIES:
            raise ValueError("Interactant A must be a moving entity.")
        
    def interactants(self):
        return f"{self.interactant_a_kind.name.capitalize()} {self.interactant_b_kind.name.capitalize()}"
        
@dataclass
class Interaction():
    interactant_a: tuple[int, int, int]
    interactant_b: tuple[int, int, int]
    interaction_kind: InteractionKind

# IMPLEMENT NEW INTERACTIONS HERE (interactant_a is automatically discarded)
def interaction_vacancy_vacancy(interactant_a, interactant_b, entities, next_entities):  # metrics could be added here as well
    next_entities.divacancy.add(choice(points_on_line(*interactant_a, *interactant_b)))
    entities.vacancy.discard(interactant_b) #Vendy's comment for Eda, shouldn't here be also the interactant_a discarded?, same for interstitials

def interaction_vacancy_divacancy(interactant_a, interactant_b, entities, next_entities):
    next_entities.trivacancy.add(interactant_b)
    next_entities.divacancy.discard(interactant_b)
    entities.divacancy.discard(interactant_b)

def interaction_vacancy_oxygen(interactant_a, interactant_b, entities, next_entities):
    next_entities.vac_oxygen.add(interactant_b)
    next_entities.oxygen.discard(interactant_b)
    entities.oxygen.discard(interactant_b)

def interaction_vacancy_vac_oxygen(interactant_a, interactant_b, entities, next_entities):
    next_entities.divac_oxygen.add(interactant_b)
    next_entities.vac_oxygen.discard(interactant_b)
    entities.vac_oxygen.discard(interactant_b)

def interaction_vacancy_interstitial(interactant_a, interactant_b, entities, next_entities):
    entities.interstitial.discard(interactant_b)

def interaction_interstitial_divacancy(interactant_a, interactant_b, entities, next_entities):
    entities.divacancy.discard(interactant_b)
    next_entities.divacancy.discard(interactant_b)
    next_entities.vacancy.add(interactant_b)

def interaction_interstitial_trivacancy(interactant_a, interactant_b, entities, next_entities):
    entities.trivacancy.discard(interactant_b)
    next_entities.trivacancy.discard(interactant_b)
    next_entities.divacancy.add(interactant_b)

def interaction_vacancy_divacancy_oxygen(interactant_a, interactant_b, entities, next_entities):
    next_entities.trivac_oxygen.add(interactant_b)
    next_entities.divac_oxygen.discard(interactant_b)
    entities.divac_oxygen.discard(interactant_b)

def interaction_interstitial_interstitial(interactant_a, interactant_b, entities, next_entities):
    next_entities.diinterstitial.add(choice(points_on_line(*interactant_a, *interactant_b)))
    entities.interstitial.discard(interactant_b)

def interaction_vacancy_diinterstitial(interactant_a, interactant_b, entities, next_entities):
    next_entities.interstitial.add(interactant_b)
    next_entities.diinterstitial.discard(interactant_b)
    entities.diinterstitial.discard(interactant_b)

def interaction_interstitial_carbon(interactant_a, interactant_b, entities, next_entities):
    next_entities.carbon_interstitial.add(interactant_b)
    next_entities.carbon.discard(interactant_b)
    entities.carbon.discard(interactant_b)

def interaction_interstitial_divacancy_oxygen(interactant_a, interactant_b, entities, next_entities):
    next_entities.vac_oxygen.add(interactant_b)
    next_entities.divac_oxygen.discard(interactant_b)
    entities.divac_oxygen.discard(interactant_b)
    
def interaction_interstitial_trivacancy_oxygen(interactant_a, interactant_b, entities, next_entities):
    next_entities.divac_oxygen.add(interactant_b)
    next_entities.trivac_oxygen.discard(interactant_b)
    entities.trivac_oxygen.discard(interactant_b)

def interaction_vacancy_carbon_interstitial(interactant_a, interactant_b, entities, next_entities):
    next_entities.carbon.add(interactant_b)
    next_entities.carbon_interstitial.discard(interactant_b)
    entities.carbon_interstitial.discard(interactant_b)

def interaction_interstitial_vacancy_oxygen(interactant_a, interactant_b, entities, next_entities):
    next_entities.oxygen.add(interactant_b)
    next_entities.vac_oxygen.discard(interactant_b)
    entities.vac_oxygen.discard(interactant_b)

INTERACTION_KINDS = [
    InteractionKind(EntityKind.VACANCY, EntityKind.VACANCY, 7.7/GRID_UNIT/10, 0.107, interaction_vacancy_vacancy),
    InteractionKind(EntityKind.VACANCY, EntityKind.DIVACANCY, 9.9/GRID_UNIT/10, 0.226, interaction_vacancy_divacancy),
    InteractionKind(EntityKind.VACANCY, EntityKind.OXYGEN, 5.0/GRID_UNIT/10, 0.029, interaction_vacancy_oxygen),
    InteractionKind(EntityKind.VACANCY, EntityKind.VAC_OXYGEN, 8.4/GRID_UNIT/10, 0.139, interaction_vacancy_vac_oxygen),
    InteractionKind(EntityKind.VACANCY, EntityKind.INTERSTITIAL, 16.0/GRID_UNIT/10, 0.956, interaction_vacancy_interstitial),
    InteractionKind(EntityKind.INTERSTITIAL, EntityKind.DIVACANCY, 15.8/GRID_UNIT/10, 0.934, interaction_interstitial_divacancy),
    InteractionKind(EntityKind.INTERSTITIAL, EntityKind.TRIVACANCY, 12.4/GRID_UNIT/10, 0.445, interaction_interstitial_trivacancy),
    InteractionKind(EntityKind.VACANCY, EntityKind.DIVAC_OXYGEN, 5.7/GRID_UNIT/10, 0.043, interaction_vacancy_divacancy_oxygen),
    InteractionKind(EntityKind.INTERSTITIAL, EntityKind.INTERSTITIAL, 7.9/GRID_UNIT/10, 0.118, interaction_interstitial_interstitial),
    InteractionKind(EntityKind.VACANCY, EntityKind.DIINTERSTITIAL, 15.3/GRID_UNIT/10, 0.849, interaction_vacancy_diinterstitial),
    InteractionKind(EntityKind.INTERSTITIAL, EntityKind.CARBON, 14.2/GRID_UNIT/10, 0.673, interaction_interstitial_carbon),
    InteractionKind(EntityKind.INTERSTITIAL, EntityKind.DIVAC_OXYGEN, 5.1/GRID_UNIT/10, 0.031, interaction_interstitial_divacancy_oxygen),
    InteractionKind(EntityKind.INTERSTITIAL, EntityKind.TRIVAC_OXYGEN, 11.7/GRID_UNIT/10, 0.374, interaction_interstitial_trivacancy_oxygen),
    InteractionKind(EntityKind.VACANCY, EntityKind.CARBON_INTERSTITIAL, 8.6/GRID_UNIT/10, 0.149, interaction_vacancy_carbon_interstitial),
    InteractionKind(EntityKind.INTERSTITIAL, EntityKind.VAC_OXYGEN, 8.6/GRID_UNIT/10, 0.149, interaction_interstitial_vacancy_oxygen),
]
# /IMPLEMENT NEW INTERACTIONS HERE

@dataclass
class PlotFormat():
    color: str
    marker: str
    label: str

# ADD FORMAT FOR NEW ENTITIES HERE
plot_format = {
    EntityKind.VACANCY: PlotFormat("cornflowerblue", "o", "$Vacancy_{init}$"),
    EntityKind.INTERSTITIAL: PlotFormat("gold", "D", "$Interstitials_{init}$"),
    EntityKind.OXYGEN: PlotFormat("paleturquoise", ".", "Oxygens"),
    EntityKind.VAC_OXYGEN: PlotFormat("turquoise", "^", "VO"),
    EntityKind.DIVAC_OXYGEN: PlotFormat("darkturquoise", "v", "V$_2$O"),
    EntityKind.DIVACANCY: PlotFormat("royalblue", "x", "V$_2$"),
    EntityKind.TRIVACANCY: PlotFormat("navy", "X", "V$_3$"),
    EntityKind.TRIVAC_OXYGEN: PlotFormat("darkblue", "P", "V$_3$O"),
    EntityKind.DIINTERSTITIAL: PlotFormat("darkgoldenrod", "d", "I$_2$"),
    EntityKind.CARBON: PlotFormat("black", "s", "Carbons"),
    EntityKind.CARBON_INTERSTITIAL: PlotFormat("gray", "P", "C$_I$"),
}

def plot_lattice(entities, step):
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    for entity_kind in EntityKind:
        if entity_kind == EntityKind.OXYGEN or entity_kind == EntityKind.CARBON:
            continue  # Skip oxygen
        if (entities_to_plot := entities.get_kind(entity_kind)):
            format = plot_format[entity_kind]
            dx, dy, dz = zip(*entities_to_plot)
            ax.scatter(dx, dy, dz, c=format.color, marker=format.marker, label=f"{len(entities_to_plot)} {format.label}")    

    # Set the limits for x, y, z axes
    ax.set_xlim(0, LATTICE_SIZE - 1)
    ax.set_ylim(0, LATTICE_SIZE - 1)
    ax.set_zlim(0, LATTICE_SIZE - 1)
    
    # Set tick labels to show only whole numbers
    ticks = range(0, LATTICE_SIZE, 10)  # Change the step if necessary for very large lattices
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_zticks(ticks)

    # Axis labels and legend
    ax.set_xlabel('X [2.35 A]')
    ax.set_ylabel('Y [2.35 A]')
    ax.set_zlabel('Z [2.35 A]')
    ax.legend()

    # Title and saving the plot
    plt.title(f'Step {step}')
    plt.savefig(f'{OUTPUT_DIR}/step_{step}.png')
    plt.close()

def calculate_neighbors(position, max_distance):
    x, y, z = position
    neighbors = []
    for dx, dy, dz in product(range(-floor(max_distance), ceil(max_distance) + 1), repeat=3):  # Range based on max_distance
        if dx == 0 and dy == 0 and dz == 0:
            continue  # Skip the current position itself
        if dx**2 + dy**2 + dz**2 <= max_distance**2:
            new_position = (x + dx, y + dy, z + dz)
            if 0 <= new_position[0] < LATTICE_SIZE and 0 <= new_position[1] < LATTICE_SIZE and 0 <= new_position[2] < LATTICE_SIZE:
                neighbors.append(new_position)
    return neighbors

def points_on_line(x0, y0, z0, x1, y1, z1):
    """
    Generate points on a line from (x0, y0, z0) to (x1, y1, z1) using a 3D Bresenham's line algorithm
    implementations from ChatGPT
    """
    points = []
    dx = abs(x1 - x0)
    dy = abs(y1 - y0)
    dz = abs(z1 - z0)
    x_sign = 1 if x0 < x1 else -1
    y_sign = 1 if y0 < y1 else -1
    z_sign = 1 if z0 < z1 else -1

    if dx >= dy and dx >= dz:        
        err_1 = 2 * dy - dx
        err_2 = 2 * dz - dx
        while x0 != x1:
            points.append((x0, y0, z0))
            if err_1 > 0:
                y0 += y_sign
                err_1 -= 2 * dx
            if err_2 > 0:
                z0 += z_sign
                err_2 -= 2 * dx
            err_1 += 2 * dy
            err_2 += 2 * dz
            x0 += x_sign
    elif dy >= dx and dy >= dz:
        err_1 = 2 * dx - dy
        err_2 = 2 * dz - dy
        while y0 != y1:
            points.append((x0, y0, z0))
            if err_1 > 0:
                x0 += x_sign
                err_1 -= 2 * dy
            if err_2 > 0:
                z0 += z_sign
                err_2 -= 2 * dy
            err_1 += 2 * dx
            err_2 += 2 * dz
            y0 += y_sign
    else:
        err_1 = 2 * dy - dz
        err_2 = 2 * dx - dz
        while z0 != z1:
            points.append((x0, y0, z0))
            if err_1 > 0:
                y0 += y_sign
                err_1 -= 2 * dz
            if err_2 > 0:
                x0 += x_sign
                err_2 -= 2 * dz
            err_1 += 2 * dy
            err_2 += 2 * dx
            z0 += z_sign

    points.append((x1, y1, z1))  # Ensure the last point is included
    return points


# Preprocess interactions
interaction_kinds_by_a_kind = {
    entity_kind: [interaction for interaction in INTERACTION_KINDS if interaction.interactant_a_kind is entity_kind]
    for entity_kind in MOVING_ENTITIES
}

# Create directory for output figures if it does not exist
OUTPUT_DIR = './DemoFigEntities'
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# Initialize entities
entities = Entities()
# Initialize displacements
num_displacements = 500
while len(entities.vacancy) < num_displacements:
    entities.vacancy.add(tuple(np.random.randint(0, LATTICE_SIZE, size=3)))

# Initialize interstitials
num_interstitials = 500
while len(entities.interstitial) < num_interstitials:
    entities.interstitial.add(tuple(np.random.randint(0, LATTICE_SIZE, size=3)))

# Initialize oxygen
SAMPLE_SIZE_X = 2.6e6  # nm
SAMPLE_SIZE_Y = 2.6e6  # nm
SAMPLE_SIZE_Z = 150e3  # nm
SAMPLE_OXYGEN_COUNT = 1e18
SAMPLE_CARBON_COUNT = 3e16

number_of_oxygens = floor(SAMPLE_OXYGEN_COUNT / (SAMPLE_SIZE_X * SAMPLE_SIZE_Y * SAMPLE_SIZE_Z) * (LATTICE_SIZE * GRID_UNIT)**3)
while len(entities.oxygen) < number_of_oxygens:
    entities.oxygen.add(tuple(np.random.randint(0, LATTICE_SIZE, size=3)))
number_of_carbons = floor(SAMPLE_CARBON_COUNT / (SAMPLE_SIZE_X * SAMPLE_SIZE_Y * SAMPLE_SIZE_Z) * (LATTICE_SIZE * GRID_UNIT)**3)
while len(entities.carbon) < number_of_carbons:
    entities.carbon.add(tuple(np.random.randint(0, LATTICE_SIZE, size=3)))
# Initialize metrics
interaction_counts = defaultdict(int)

class EntityCounts(EntitiesMeta):
    container = list

    def injest(self, **kwargs):
        for key, value in kwargs.items():
            getattr(self, key).append(value)

entity_counts = EntityCounts()
entity_counts.injest(**entities.report_counts())

for step in tqdm(range(STEPS)):
    plot_lattice(entities, step)

    next_entities = entities.copy()
    for entity_kind in MOVING_ENTITIES:
        next_entities.get_kind(entity_kind).clear()

    for entity_kind in MOVING_ENTITIES:
        entities_to_interact = entities.get_kind(entity_kind)
        while entities_to_interact:
            interactant_a = entities_to_interact.pop()

            interactions = []
            for interaction_kind in interaction_kinds_by_a_kind[entity_kind]:
                neighbouring_interactant_b = [
                    entity for entity in calculate_neighbors(interactant_a, interaction_kind.radius) 
                    if entity in entities.get_kind(interaction_kind.interactant_b_kind)
                ]
                interactions.extend(Interaction(interactant_a, interactant_b, interaction_kind) for interactant_b in neighbouring_interactant_b)

            shuffle(interactions)

            for interaction in interactions:
                if random() < interaction.interaction_kind.probability:
                    interaction.interaction_kind.process(
                        interaction.interactant_a,
                        interaction.interactant_b,
                        entities,
                        next_entities,
                    )
                    # log interaction counts
                    interaction_counts[interaction.interaction_kind.interactants()] += 1
                    break
            else:
                next_entities.get_kind(entity_kind).add(interactant_a)


    entities = next_entities

    # Perform random walks
    for entity_kind in MOVING_ENTITIES:
        entities_to_move = entities.get_kind(entity_kind)
        moved_entities = set()
        for entity in entities_to_move:
            move = np.random.randint(-1, 2, size=3)
            new_position = tuple(np.array(entity) + move)
            if 0 <= min(new_position) < LATTICE_SIZE and max(new_position) < LATTICE_SIZE:
                if new_position not in moved_entities:  # removed checks for presence of other entities in the same coordinates
                    moved_entities.add(new_position)
                else:
                    moved_entities.add(entity)
            # Boundary or invalid move cases are ignored (displacement lost)
        entities_to_move.clear()
        entities_to_move.update(moved_entities)

    # Record counts
    entity_counts.injest(**entities.report_counts())
    
plot_lattice(entities, STEPS)

# Final counts plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

for entity_kind in EntityKind:
    if entity_kind == EntityKind.OXYGEN or entity_kind == EntityKind.CARBON:
        continue
    format = plot_format[entity_kind]
    ax1.plot(range(STEPS + 1), entity_counts.get_kind(entity_kind), label=format.label)
ax1.set_xlabel('Step')
ax1.set_ylabel('Count')
ax1.set_title('Simulation Results Over Time')
ax1.legend()
print(interaction_counts.keys())
key_mapping = {
    'Vacancy Vacancy': 'V+V',
    'Vacancy Interstitial': 'V+I',
    'Interstitial Divacancy': 'I+V2',
    'Vacancy Divacancy': 'V+V2',
    'Vacancy Oxygen': 'V+O',
    'Vacancy Vac_oxygen': 'V+VO',
    'Vacancy Divac_oxygen': 'V+V2O',
    'Vacancy Diinterstitial': 'V+I2',
    'Interstitial Interstitial': 'I+I',
    'Interstitial Divacancy': 'I+V2',
    'Interstitial Trivacancy': 'I+V3',
    'Interstitial_Divac_oxygen': 'I+V2O',
    'Interstitial Carbon': 'C+I'

    
}

# Update the interaction_counts dictionary to use the new abbreviated keys
abbreviated_interaction_counts = {key_mapping[key]: value for key, value in interaction_counts.items() if key in key_mapping}

# Now plot using the updated dictionary
ax2.bar(abbreviated_interaction_counts.keys(), abbreviated_interaction_counts.values(), color='tab:blue')
ax2.set_title('Interaction Counts')
ax2.set_title('Interaction counts')

plt.savefig(f'{OUTPUT_DIR}/final_counts.png')
plt.tight_layout()

plt.show()
