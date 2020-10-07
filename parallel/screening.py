from pymatgen.io.cif import CifParser
import candidates
import two_bond
import three_bond
import four_bond

def screening(cif_file, neighbor_distance, lengths, length, arm_length, upper_limit, lower_limit, step_size, angle_limit, mof, n_mofs):
    parser = CifParser(cif_file)
    structure = parser.get_structures(primitive=False)[0]
    screening_candidates = candidates.candidates(structure, neighbor_distance)
    screening_candidates = candidates.ring_checker(structure, screening_candidates, neighbor_distance)

    if len(screening_candidates) == 0:
        return ""

    if len(lengths) == 1:
        return two_bond.two_bond_screening(cif_file, structure, screening_candidates, length, upper_limit,
                            lower_limit, angle_limit, step_size, mof, n_mofs)

    if len(lengths) == 3:
        return three_bond.three_bond_screening(cif_file, structure, screening_candidates, lengths, arm_length, upper_limit,
                            lower_limit, angle_limit, step_size, mof, n_mofs)
    elif len(lengths) == 4:
        return four_bond.four_bond_screening(cif_file, structure, screening_candidates, lengths, arm_length, upper_limit,
                            lower_limit, angle_limit, step_size, mof, n_mofs)
    else:
        return ""
