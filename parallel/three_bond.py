import filters
import numpy as np
import ntpath

# THREE BOND SCREENING
def three_bond_screening(filename, structure, candidates, lengths, arm_length, upper_limit, lower_limit, angle_limit, step_size, mof, n_mofs):
    output = []
    #print(file_path)
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

                #print("3-Bond Screening, checking MOF {0}/{1}: MOF {2}, site {3}, {4} and {5}".format(mof, n_mofs, filename, candidates[i][2], candidates[j][2], candidates[k][2]), end='\r')

                center = (h_1_pos + h_2_pos + h_3_pos)/3

                initial_distance_1 = filters.distance(h_1_pos, center)
                initial_distance_2 = filters.distance(h_2_pos, center)
                initial_distance_3 = filters.distance(h_3_pos, center)

                h_distances = [initial_distance_1, initial_distance_2, initial_distance_3]

                if filters.common_member(candidates[i], candidates[j], candidates[k]):
                    continue

                if any(h_distance <= arm_length*0.65 for h_distance in h_distances):
                    continue

                if any(h_distance >= arm_length*1.5 for h_distance in h_distances):
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

                    h_1_pos = filters.rotate(h_1_pos, axis_1, shift_1, step_size_1)
                    c_1_pos = filters.rotate(h_1_pos, axis_1, shift_1, step_size_1)
                    temp_distance = filters.distance(center, h_1_pos)

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

                    h_2_pos = filters.rotate(h_2_pos, axis_2, shift_2, step_size_2)
                    c_2_pos = filters.rotate(h_2_pos, axis_2, shift_2, step_size_2)
                    temp_distance = filters.distance(center, h_2_pos)

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

                    h_3_pos = filters.rotate(h_3_pos, axis_3, shift_3, step_size_3)
                    c_3_pos = filters.rotate(h_3_pos, axis_3, shift_3, step_size_3)
                    temp_distance = filters.distance(center, h_3_pos)

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
                angle1 = filters.angle(hc_vector_1, cc_vector_1)

                cc_vector_2 = center - c_2_pos
                hc_vector_2 = h_2_pos - c_2_pos
                angle2 = filters.angle(hc_vector_2, cc_vector_2)

                cc_vector_3 = center - c_3_pos
                hc_vector_3 = h_3_pos - c_3_pos
                angle3 = filters.angle(hc_vector_3, cc_vector_3)

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

                output.append("HIT: Site {0}, Site {1}, Site {2} and Site {3} of MOF {4} allow cross-linking length of {5}!".format(candidates[i][0], candidates[j][0], candidates[k][0], ntpath.basename(filename), arm_length))

    return output
