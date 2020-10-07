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
