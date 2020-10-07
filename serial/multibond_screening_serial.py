import os
import numpy as np
import time
import datetime
from pymatgen.io.cif import CifParser


database_path = "path"
length = 0.
arm_length = 6.528  # only if cross-linker is symmetrical
lengths = np.array([arm_length, arm_length, arm_length])
upper_limit = 1.7
lower_limit = 1.6
step_size = 1000
angle_limit = 180.0
neighbor_distance = 1.55
n_mofs = len([name for name in os.listdir(database_path) if name.endswith(".cif")])
mof = 0
startTime = time.time()
output = "core_screening_" + datetime.datetime.now().strftime("%I:%M:%S-%d%m%Y") + ".log"


# MODULES FOR CALCULATION
def angle(v1, v2):
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    return 180*(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))/np.pi


def distance(vector_1, vector_2):
    return np.linalg.norm(vector_2 - vector_1)


def rotation(vector, axis, theta):
    return vector*np.cos(-theta)+np.cross(axis, vector)*np.sin(-theta)+axis*np.dot(axis, vector)*(1-np.cos(-theta))


def rotate(point, axis, shift, step_size):
    theta = np.pi/step_size
    point = point - shift
    point = rotation(point, axis, theta)
    point = point + shift

    return point


def common_member(candidate_1, candidate_2, candidate_3, candidate_4 = (0, 0, 0, 0)):
    set_1 = set(candidate_1)
    set_2 = set(candidate_2)
    set_3 = set(candidate_3)
    set_4 = set(candidate_4)
    if (set_1 & set_2 or set_1 & set_3 or set_1 & set_4 or set_2 & set_3 or set_2 & set_4 or set_3 & set_4):
        return True
    else:
        return False


def symmetry(liste):
    average = np.average(liste)

    if any(distance <= average*0.9 or distance >= average*1.1 for distance in liste):
        return False
    else:
        return True

# GET CANDIDATE
def candidates(structure, radius):
    hydrogens = []
    candidates = []

    for ndx, atom in enumerate(structure):
        if atom.species_string != "H":
            continue
        else:
            hydrogens.append(ndx)

    for hydrogen in hydrogens:
        neighbors = structure.get_neighbors(structure[hydrogen], 1.5, True, False)
        if structure[neighbors[0][2]].species_string == "C" and len(neighbors) == 1:
            candidates.append((hydrogen, neighbors[0][2]))
    return candidates


# FILTER CANDIDATES
def ring_checker(structure, candidates, neighbor_distance):
    filtered_candidates = []
    #print("Im Ring-Checker 1: ", len(candidates))
    for hydrogens, carbons in candidates:
        neighbor_list = []
        next_neighbor_list_1 = []
        next_neighbor_list_2 = []
        temp_neighbors = structure.get_neighbors(structure[carbons], neighbor_distance, True, False)

        #print("Nachbarliste der Kandidaten: ", len(temp_neighbors), z)

        if len(temp_neighbors) != 3:
            #print("continue")
            continue
        else:
            neighbor_list.append((temp_neighbors[0][0].species_string, temp_neighbors[0][2]))
            neighbor_list.append((temp_neighbors[1][0].species_string, temp_neighbors[1][2]))
            neighbor_list.append((temp_neighbors[2][0].species_string, temp_neighbors[2][2]))

        temp_neighbor_list = neighbor_list

        #print(temp_neighbor_list)

        for element in neighbor_list:
            #print(element)
            if element[0] == 'H':
                temp_neighbor_list.remove(element)

        neighbor_list = temp_neighbor_list
        del temp_neighbor_list

        #print(neighbor_list, z)

        if len(neighbor_list) != 2:
            continue
        else:
            #print("pass bois")
            pass

        if neighbor_list[0][0] != 'C' or neighbor_list[1][0] != 'C':
            #print("Fremdatome")
            continue
        else:
            pass

        direction_1 = structure.get_neighbors(structure[neighbor_list[0][1]], neighbor_distance, True, False)
        direction_2 = structure.get_neighbors(structure[neighbor_list[1][1]], neighbor_distance, True, False)

        #print(direction_1, direction_2)

        if len(direction_1) != 3 or len(direction_2) != 3:
            #print(len(direction_1), len(direction_2))
            continue
        else:
            next_neighbor_list_1.append((direction_1[0][0].species_string, direction_1[0][2]))
            next_neighbor_list_1.append((direction_1[1][0].species_string, direction_1[1][2]))
            next_neighbor_list_1.append((direction_1[2][0].species_string, direction_1[2][2]))
            next_neighbor_list_2.append((direction_2[0][0].species_string, direction_2[0][2]))
            next_neighbor_list_2.append((direction_2[1][0].species_string, direction_2[1][2]))
            next_neighbor_list_2.append((direction_2[2][0].species_string, direction_2[2][2]))

        #ERSTER ACHSENPUNKT
        c_counter = 0
        for element in next_neighbor_list_1:
            if element[0] == 'C':
                c_counter = c_counter + 1

        if c_counter == 3:
            axis_1 = neighbor_list[0][1]
        else:
            next_next_neighbor_list = []
            for element in next_neighbor_list_1:
                if element[0] == 'H':
                    next_neighbor_list_1.remove(element)
            for element in next_neighbor_list_1:
                if element[1] == carbons:
                    next_neighbor_list_1.remove(element)

            if len(next_neighbor_list_1) != 1:
                continue
            else:
                pass
            if next_neighbor_list_1[0][0] != 'C':
                continue
            else:
                pass

            next_direction = structure.get_neighbors(structure[next_neighbor_list_1[0][1]], neighbor_distance, True, False)

            if len(next_direction) != 3:
                continue
            else:
                pass

            next_next_neighbor_list.append((next_direction[0][0].species_string, next_direction[0][2]))
            next_next_neighbor_list.append((next_direction[1][0].species_string, next_direction[1][2]))
            next_next_neighbor_list.append((next_direction[2][0].species_string, next_direction[2][2]))

            c_counter = 0

            for element in next_next_neighbor_list:
                if element[0] == 'C':
                    c_counter = c_counter + 1

            if c_counter != 3:
                continue
            else:
                axis_1 = next_neighbor_list_1[0][1]

        #ZWEITER ACHSENPUNKT
        c_counter = 0
        for element in next_neighbor_list_2:
            if element[0] == 'C':
                c_counter = c_counter + 1

        if c_counter == 3:
            axis_2 = neighbor_list[1][1]
        else:
            next_next_neighbor_list = []
            for element in next_neighbor_list_2:
                if element[0] == 'H':
                    next_neighbor_list_2.remove(element)
            for element in next_neighbor_list_2:
                if element[1] == carbons:
                    next_neighbor_list_2.remove(element)

            if len(next_neighbor_list_2) != 1:
                continue
            else:
                pass

            if next_neighbor_list_2[0][0] != 'C':
                continue
            else:
                pass

            next_direction = structure.get_neighbors(structure[next_neighbor_list_2[0][1]], neighbor_distance, True, False)

            if len(next_direction) != 3:
                continue
            else:
                pass

            next_next_neighbor_list.append((next_direction[0][0].species_string, next_direction[0][2]))
            next_next_neighbor_list.append((next_direction[1][0].species_string, next_direction[1][2]))
            next_next_neighbor_list.append((next_direction[2][0].species_string, next_direction[2][2]))

            c_counter = 0

            for element in next_next_neighbor_list:
                if element[0] == 'C':
                    c_counter = c_counter + 1

            if c_counter != 3:
                continue
            else:
                axis_2 = next_neighbor_list_1[0][1]

        filtered_candidates.append((hydrogens, carbons, axis_1, axis_2))
    return filtered_candidates


# THREE BOND SCREENING
def three_bond_screening(filename, structure, candidates, lengths, upper_limit, lower_limit, angle_limit, step_size, ouput, mof, n_mofs):
    for i in range(len(candidates)):
        axis_1 = structure[candidates[i][2]].coords - structure[candidates[i][3]].coords
        axis_1 = [round(elem, 8) for elem in axis_1]
        axis_1 = np.asarray(axis_1, dtype=np.float32)
        axis_1 = axis_1 / np.linalg.norm(axis_1)

        shift_1 = structure[candidates[i][2]].coords
        shift_1 = np.asarray(shift_1, dtype=np.float32)
        h_1_pos = structure[candidates[i][0]].coords
        h_1_pos = np.asarray(h_1_pos, dtype=np.float32)
        c_1_pos = structure[candidates[i][1]].coords
        c_1_pos = np.asarray(c_1_pos, dtype=np.float32)

        for j in range(i):
            axis_2 = structure[candidates[j][2]].coords - structure[candidates[j][3]].coords
            axis_2 = [round(elem, 8) for elem in axis_2]
            axis_2 = np.asarray(axis_2, dtype=np.float32)
            axis_2 = axis_2 / np.linalg.norm(axis_2)

            shift_2 = structure[candidates[j][2]].coords
            shift_2 = np.asarray(shift_2, dtype=np.float32)
            h_2_pos = structure[candidates[j][0]].coords
            h_2_pos = np.asarray(h_2_pos, dtype=np.float32)
            c_2_pos = structure[candidates[j][1]].coords
            c_2_pos = np.asarray(c_2_pos, dtype=np.float32)

            for k in range(j):
                axis_3 = structure[candidates[k][2]].coords - structure[candidates[k][3]].coords
                axis_3 = [round(elem, 8) for elem in axis_3]
                axis_3 = np.asarray(axis_3, dtype=np.float32)
                axis_3 = axis_3 / np.linalg.norm(axis_3)

                shift_3 = structure[candidates[k][2]].coords
                shift_3 = np.asarray(shift_3, dtype=np.float32)
                h_3_pos = structure[candidates[k][0]].coords
                h_3_pos = np.asarray(h_3_pos, dtype=np.float32)
                c_3_pos = structure[candidates[k][1]].coords
                c_3_pos = np.asarray(c_3_pos, dtype=np.float32)

                print("3-Bond Screening, checking MOF {0}/{1}: MOF {2}, site {3}, {4} and {5}".format(mof, n_mofs, filename, candidates[i][2], candidates[j][2], candidates[k][2]), end='\r')

                center = (h_1_pos + h_2_pos + h_3_pos)/3

                initial_distance_1 = distance(h_1_pos, center)
                initial_distance_2 = distance(h_2_pos, center)
                initial_distance_3 = distance(h_3_pos, center)

                h_distances = [initial_distance_1, initial_distance_2, initial_distance_3]

                if common_member(candidates[i], candidates[j], candidates[k]):
                    continue

                if any(h_distance <= arm_length*0.3 for h_distance in h_distances):
                    continue

                counter_1 = 0
                counter_2 = 0
                counter_3 = 0

                step_size_1 = step_size
                step_size_2 = step_size
                step_size_3 = step_size

                h_distance = initial_distance_1
                temp_distance = 0.

                while(h_distance >= temp_distance):
                    #print("Rotating Linker 1", end='\r')

                    h_1_pos = rotate(h_1_pos, axis_1, shift_1, step_size_1)
                    c_1_pos = rotate(h_1_pos, axis_1, shift_1, step_size_1)
                    temp_distance = distance(center, h_1_pos)

                    if counter_1 == 0:
                        if temp_distance > h_distance:
                            step_size_1 = -step_size
                        else:
                            pass
                    else:
                        pass

                    if temp_distance > h_distance:
                        break

                    h_distance = temp_distance
                    counter_1 += 1

                final_distance_1 = h_distance
                h_distance = initial_distance_2
                temp_distance = 0.

                while(h_distance >= temp_distance):
                    #print("Rotating Linker 2", end='\r')

                    h_2_pos = rotate(h_2_pos, axis_2, shift_2, step_size_2)
                    c_2_pos = rotate(h_2_pos, axis_2, shift_2, step_size_2)
                    temp_distance = distance(center, h_2_pos)

                    if counter_2 == 0:
                        if temp_distance > h_distance:
                            step_size_2 = -step_size
                        else:
                            pass
                    else:
                        pass

                    if temp_distance > h_distance:
                        break

                    h_distance = temp_distance
                    counter_2 += 1

                final_distance_2 = h_distance
                h_distance = initial_distance_3
                temp_distance = 0.

                while(h_distance >= temp_distance):
                    #print("Rotating Linker 3", end='\r')

                    h_3_pos = rotate(h_3_pos, axis_3, shift_3, step_size_3)
                    c_3_pos = rotate(h_3_pos, axis_3, shift_3, step_size_3)
                    temp_distance = distance(center, h_3_pos)

                    if counter_3 == 0:
                        if temp_distance > h_distance:
                            step_size_3 = -step_size
                        else:
                            pass
                    else:
                        pass

                    if temp_distance > h_distance:
                        break

                    h_distance = temp_distance
                    counter_3 += 1

                final_distance_3 = h_distance

                new_center = (h_1_pos + h_2_pos + h_3_pos)/3

                if(distance(center, new_center) >= 0.8*arm_length):
                    continue

                distances = np.array([final_distance_1, final_distance_2, final_distance_3])

                sorted_l = np.sort(lengths)
                sorted_d = np.sort(distances)

                if symmetry(sorted_d):
                    pass
                else:
                    continue

                cc_vector_1 = center - c_1_pos
                hc_vector_1 = h_1_pos - c_1_pos
                angle1 = angle(hc_vector_1, cc_vector_1)

                cc_vector_2 = center - c_2_pos
                hc_vector_2 = h_2_pos - c_2_pos
                angle2 = angle(hc_vector_2, cc_vector_2)

                cc_vector_3 = center - c_3_pos
                hc_vector_3 = h_3_pos - c_3_pos
                angle3 = angle(hc_vector_3, cc_vector_3)

                if angle1 >= angle_limit:
                    continue
                elif angle2 >= angle_limit:
                    continue
                elif angle3 >= angle_limit:
                    continue
                else:
                    pass

                if np.average(sorted_d) <= np.average(sorted_l)*upper_limit and np.average(sorted_d) >= np.average(sorted_l)*lower_limit:
                        pass
                else:
                    continue

                print("HIT: Site {0}, Site {1} and Site {2} of MOF {3} allow cross-linking length of {4}!\n".format(candidates[i][0], candidates[j][0], candidates[k][0], filename, arm_length))
                with open(output, 'a+') as the_file:
                    the_file.write("HIT: Site {0}, {1} and {2} of MOF {3} allow cross-linking lengths of {4}!\n".format(candidates[i][0], candidates[j][0], candidates[k][0], filename, arm_length))


# FOUR BOND SCREENING
def four_bond_screening(filename, structure, candidates, lengths, upper_limit, lower_limit, angle_limit, step_size, ouput, mof, n_mofs):
    for i in range(len(candidates)):
        axis_1 = structure[candidates[i][2]].coords - structure[candidates[i][3]].coords
        axis_1 = [round(elem, 8) for elem in axis_1]
        axis_1 = np.asarray(axis_1, dtype=np.float32)
        axis_1 = axis_1 / np.linalg.norm(axis_1)

        shift_1 = structure[candidates[i][2]].coords
        shift_1 = np.asarray(shift_1, dtype=np.float32)
        h_1_pos = structure[candidates[i][0]].coords
        h_1_pos = np.asarray(h_1_pos, dtype=np.float32)
        c_1_pos = structure[candidates[i][1]].coords
        c_1_pos = np.asarray(c_1_pos, dtype=np.float32)

        for j in range(i):
            axis_2 = structure[candidates[j][2]].coords - structure[candidates[j][3]].coords
            axis_2 = [round(elem, 8) for elem in axis_2]
            axis_2 = np.asarray(axis_2, dtype=np.float32)
            axis_2 = axis_2 / np.linalg.norm(axis_2)

            shift_2 = structure[candidates[j][2]].coords
            shift_2 = np.asarray(shift_2, dtype=np.float32)
            h_2_pos = structure[candidates[j][0]].coords
            h_2_pos = np.asarray(h_2_pos, dtype=np.float32)
            c_2_pos = structure[candidates[j][1]].coords
            c_2_pos = np.asarray(c_2_pos, dtype=np.float32)

            for k in range(j):
                axis_3 = structure[candidates[k][2]].coords - structure[candidates[k][3]].coords
                axis_3 = [round(elem, 8) for elem in axis_3]
                axis_3 = np.asarray(axis_3, dtype=np.float32)
                axis_3 = axis_3 / np.linalg.norm(axis_3)

                shift_3 = structure[candidates[k][2]].coords
                shift_3 = np.asarray(shift_3, dtype=np.float32)
                h_3_pos = structure[candidates[k][0]].coords
                h_3_pos = np.asarray(h_3_pos, dtype=np.float32)
                c_3_pos = structure[candidates[k][1]].coords
                c_3_pos = np.asarray(c_3_pos, dtype=np.float32)

                for l in range(k):
                    axis_4 = structure[candidates[l][2]].coords - structure[candidates[l][3]].coords
                    axis_4 = [round(elem, 8) for elem in axis_4]
                    axis_4 = np.asarray(axis_4, dtype=np.float32)
                    axis_4 = axis_4 / np.linalg.norm(axis_4)

                    shift_4 = structure[candidates[l][2]].coords
                    shift_4 = np.asarray(shift_4, dtype=np.float32)
                    h_4_pos = structure[candidates[l][0]].coords
                    h_4_pos = np.asarray(h_4_pos, dtype=np.float32)
                    c_4_pos = structure[candidates[l][1]].coords
                    c_4_pos = np.asarray(c_4_pos, dtype=np.float32)

                    print("4-Bond Screening, checking MOF {0}/{1}: MOF {2}, site {3}, {4}, {5} and {6}".format(mof, n_mofs, filename, candidates[i][0], candidates[j][0], candidates[k][0], candidates[l][0]), end='\r')

                    center = (h_1_pos + h_2_pos + h_3_pos + h_4_pos)/4

                    initial_distance_1 = distance(h_1_pos, center)
                    initial_distance_2 = distance(h_2_pos, center)
                    initial_distance_3 = distance(h_3_pos, center)
                    initial_distance_4 = distance(h_4_pos, center)

                    if common_member(candidates[i], candidates[j], candidates[k], candidates[l]):
                        continue

                    initial_distances = [initial_distance_1, initial_distance_2, initial_distance_3, initial_distance_4]

                    if any(initial_distance <= arm_length*0.65 for initial_distance in initial_distances):
                        continue

                    counter_1 = 0
                    counter_2 = 0
                    counter_3 = 0
                    counter_4 = 0

                    step_size_1 = step_size
                    step_size_2 = step_size
                    step_size_3 = step_size
                    step_size_4 = step_size

                    h_distance = initial_distance_1
                    temp_distance = 0.

                    while(h_distance >= temp_distance):
                        #print("Rotating Linker 1", end='\r')

                        h_1_pos = rotate(h_1_pos, axis_1, shift_1, step_size_1)
                        c_1_pos = rotate(h_1_pos, axis_1, shift_1, step_size_1)
                        temp_distance = distance(center, h_1_pos)

                        if counter_1 == 0:
                            if temp_distance > h_distance:
                                step_size_1 = -step_size
                            else:
                                pass
                        else:
                            pass

                        if temp_distance > h_distance:
                            break

                        h_distance = temp_distance
                        counter_1 += 1

                    final_distance_1 = h_distance
                    h_distance = initial_distance_2
                    temp_distance = 0.

                    while(h_distance >= temp_distance):
                        #print("Rotating Linker 2", end='\r')

                        h_2_pos = rotate(h_2_pos, axis_2, shift_2, step_size_2)
                        c_2_pos = rotate(h_2_pos, axis_2, shift_2, step_size_2)
                        temp_distance = distance(center, h_2_pos)

                        if counter_2 == 0:
                            if temp_distance > h_distance:
                                step_size_2 = -step_size
                            else:
                                pass
                        else:
                            pass

                        if temp_distance > h_distance:
                            break

                        h_distance = temp_distance
                        counter_2 += 1

                    final_distance_2 = h_distance
                    h_distance = initial_distance_3
                    temp_distance = 0.

                    while(h_distance >= temp_distance):
                        #print("Rotating Linker 3", end='\r')

                        h_3_pos = rotate(h_3_pos, axis_3, shift_3, step_size_3)
                        c_3_pos = rotate(h_3_pos, axis_3, shift_3, step_size_3)
                        temp_distance = distance(center, h_3_pos)

                        if counter_3 == 0:
                            if temp_distance > h_distance:
                                step_size_3 = -step_size
                            else:
                                pass
                        else:
                            pass

                        if temp_distance > h_distance:
                            break

                        h_distance = temp_distance
                        counter_3 += 1

                    final_distance_3 = h_distance
                    h_distance = initial_distance_4
                    temp_distance = 0.

                    while(h_distance >= temp_distance):
                        #print("Rotating Linker 4", end='\r')

                        h_4_pos = rotate(h_4_pos, axis_4, shift_4, step_size_4)
                        c_4_pos = rotate(h_4_pos, axis_4, shift_4, step_size_4)
                        temp_distance = distance(center, h_4_pos)

                        if counter_4 == 0:
                            if temp_distance > h_distance:
                                step_size_4 = -step_size
                            else:
                                pass
                        else:
                            pass

                        if temp_distance > h_distance:
                            break

                        h_distance = temp_distance
                        counter_4 += 1

                    final_distance_4 = h_distance

                    new_center = (h_1_pos + h_2_pos + h_3_pos + h_4_pos)/4

                    if(distance(center, new_center) >= 0.8*arm_length):
                        continue

                    distances = np.array([final_distance_1, final_distance_2, final_distance_3, final_distance_4])

                    sorted_l = np.sort(lengths)
                    sorted_d = np.sort(distances)

                    if symmetry(sorted_d):
                        pass
                    else:
                        continue

                    cc_vector_1 = center - c_1_pos
                    hc_vector_1 = h_1_pos - c_1_pos
                    angle1 = angle(hc_vector_1, cc_vector_1)

                    cc_vector_2 = center - c_2_pos
                    hc_vector_2 = h_2_pos - c_2_pos
                    angle2 = angle(hc_vector_2, cc_vector_2)

                    cc_vector_3 = center - c_3_pos
                    hc_vector_3 = h_3_pos - c_3_pos
                    angle3 = angle(hc_vector_3, cc_vector_3)

                    cc_vector_4 = center - c_4_pos
                    hc_vector_4 = h_4_pos - c_4_pos
                    angle4 = angle(hc_vector_4, cc_vector_4)

                    if angle1 >= angle_limit:
                        continue
                    elif angle2 >= angle_limit:
                        continue
                    elif angle3 >= angle_limit:
                        continue
                    elif angle4 >= angle_limit:
                        continue
                    else:
                        pass

                    if np.average(sorted_d) <= np.average(sorted_l)*upper_limit and np.average(sorted_d) >= np.average(sorted_l)*lower_limit:
                        pass
                    else:
                        continue

                    print("HIT: Site {0}, Site {1}, Site {2} and Site {3} of MOF {4} allow cross-linking length of {5}!\n".format(candidates[i][0], candidates[j][0], candidates[k][0], candidates[l][0], filename, arm_length))
                    with open(output, 'a+') as the_file:
                        the_file.write("HIT: Site {0}, {1}, {2} and {3} of MOF {4} allow cross-linking lengths of {5}!\n".format(candidates[i][0], candidates[j][0], candidates[k][0], candidates[l][0], filename, arm_length))


if __name__ == "__main__":
    for filename in os.listdir(database_path):
        mof = mof + 1

        if filename.endswith(".cif"):

            file_path = os.path.join(database_path, filename)
            parser = CifParser(file_path)
            structure = parser.get_structures(primitive=False)[0]

            screening_candidates = candidates(structure, neighbor_distance)

            screening_candidates = ring_checker(structure, screening_candidates, neighbor_distance)
            print(len(screening_candidates))

            if len(screening_candidates) == 0:
                continue

            if len(lengths) == 3:
                three_bond_screening(filename, structure, screening_candidates, lengths, upper_limit,
                                     lower_limit, angle_limit, step_size, output, mof, n_mofs)
            elif len(lengths) == 4:
                four_bond_screening(filename, structure, screening_candidates, lengths, upper_limit,
                                    lower_limit, angle_limit, step_size, output, mof, n_mofs)
            else:
                "ERROR: WRONG LINKER TYPE!"

            with open(output, 'a') as the_file:
                the_file.write("\n")
        else:
            continue

    print("This screening took {0} seconds to execute!".format(time.time() - startTime))
