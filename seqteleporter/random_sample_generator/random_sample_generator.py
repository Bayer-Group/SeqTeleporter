import random
from typing import Tuple, List, Any, Union, Optional

# List of single-letter codes for the 20 standard amino acids
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'


def generate_random_amino_acid_sequence(min_length: int, max_length: int, amino_acids: str) -> str:
    # Ensure the length is within the specified range
    length = random.randint(min_length, max_length)

    # Generate a random sequence of amino acids
    sequence = ''.join(random.choices(amino_acids, k=length))

    return sequence


def aa_mutation_generator(aa_seq: str,
                          min_number_of_positions: int, max_number_of_positions: int,
                          min_variations_per_position: int, max_variations_per_position: int,
                          amino_acids: str) -> list:
    number_of_positions = random.randint(min_number_of_positions, max_number_of_positions)
    mutation_positions = random.sample(range(len(aa_seq)), k=number_of_positions)

    mutations = []
    for position in mutation_positions:
        wt_aa = aa_seq[position]
        number_of_variations = random.randint(min_variations_per_position, max_variations_per_position)
        variations = random.sample(amino_acids.replace(wt_aa, ""), k=number_of_variations)
        mutations.append({'position': position + 1,  # adjust to 1-indexing
                          'aa': variations})

    return sorted(mutations, key=lambda x: (x['position']))


def linked_mutation_generator(aa_seq: str, max_positions_per_set: int, max_number_of_mutation_sets: int,
                              amino_acids: str) -> Union[List[Any], None]:
    if max_number_of_mutation_sets == 0:
        return None
    number_of_mutation_sets = random.randint(1, max_number_of_mutation_sets)
    mutation_sets: list = []
    seen_positions: list = []
    for mutation_set_idx in range(0, number_of_mutation_sets):
        number_of_positions_in_this_set = random.randint(2, max_positions_per_set)
        mutation_positions = random.sample(range(len(aa_seq)), k=number_of_positions_in_this_set)
        # if any of the positions is already seen in other mutation sets, re-sample.
        while not set(mutation_positions).isdisjoint(seen_positions):
            mutation_positions = random.sample(range(len(aa_seq)), k=number_of_positions_in_this_set)
        seen_positions = seen_positions + mutation_positions
        mutation_set = []
        for position in mutation_positions:
            wt_aa = aa_seq[position]
            mut_aa = random.sample(amino_acids.replace(wt_aa, ""), k=1)[0]
            mutation_set.append((wt_aa, position, mut_aa))
        mutation_sets.append(tuple(mutation_set))

    return sorted(mutation_sets)


def random_sample_generator(min_aa_length: int = 30, max_aa_length: int = 500,
                            min_number_of_positions: int = 3, max_number_of_positions: int = 6,
                            min_variations_per_position: int = 1, max_variations_per_position: int = 10,
                            max_positions_per_linked_mutation_set: int = 3,
                            max_number_of_mutation_linked_mutation_sets: int = 3) ->\
        Tuple[str, List[Any], Optional[List[Any]]]:
    random_aa_seq = generate_random_amino_acid_sequence(min_aa_length, max_aa_length, amino_acids=AMINO_ACIDS)
    print(f"Random amino acid sequence of length {len(random_aa_seq)}: {random_aa_seq}")

    mutations = aa_mutation_generator(aa_seq=random_aa_seq,
                                      min_number_of_positions=min_number_of_positions,
                                      max_number_of_positions=max_number_of_positions,
                                      min_variations_per_position=min_variations_per_position,
                                      max_variations_per_position=max_variations_per_position,
                                      amino_acids=AMINO_ACIDS)
    positions_of_mutations = [d['position'] for d in mutations]
    print(f"\nRandom 1-indexed mutations at {len(mutations)} positions: \n {mutations}")

    linked_mutations = linked_mutation_generator(
        aa_seq=random_aa_seq,
        max_positions_per_set=max_positions_per_linked_mutation_set,
        max_number_of_mutation_sets=max_number_of_mutation_linked_mutation_sets,
        amino_acids=AMINO_ACIDS
    )

    if linked_mutations:
        positions_of_linked_mutations = [mut[1] for mut_set in linked_mutations for mut in mut_set]
        loop = 0
        while not set(positions_of_linked_mutations).isdisjoint(positions_of_mutations):
            loop += 1
            linked_mutations = linked_mutation_generator(
                aa_seq=random_aa_seq,
                max_positions_per_set=max_positions_per_linked_mutation_set,
                max_number_of_mutation_sets=max_number_of_mutation_linked_mutation_sets,
                amino_acids=AMINO_ACIDS
            )
            if loop == 100:
                linked_mutations = None
                print(f"WARNING! Unable to find linked mutations that are not overlap with mutation positions.\n"
                      f"Set linked_mutations=None")
                break
    print(f"\nRandom 1-indexed linked mutations: \n {linked_mutations}")

    return random_aa_seq, mutations, linked_mutations
