import numpy as np
import filters
import ntpath

def two_bond_screening(filename, structure, candidates, length, upper_limit, lower_limit, angle_limit, step_size, mof, n_mofs):
    output = []
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

            h_distance = filters.distance(h_2_pos, h_1_pos)
            temp_distance = 0.
            counter_1 = 0
            counter_2 = 0
            h_distance_1 = h_distance
            h_distance_2 = h_distance

            while(h_distance >= temp_distance):
                step_size_1 = step_size
                step_size_2 = step_size
                temp_distance_1 = 0.
                temp_distance_2 = 0.

                while(h_distance_1 > temp_distance_1):
                    h_1_pos = filters.rotate(h_1_pos, axis_1, shift_1, step_size_1)
                    c_1_pos = filters.rotate(c_1_pos, axis_1, shift_1, step_size_1)
                    temp_distance_1 = filters.distance(h_2_pos, h_1_pos)

                    if counter_1 == 0:
                        if temp_distance_1 > h_distance:
                            step_size_1 = -step_size
                            h_1_pos = filters.rotate(h_1_pos, axis_1, shift_1, step_size_1)
                            c_1_pos = filters.rotate(c_1_pos, axis_1, shift_1, step_size_1)
                            temp_distance_1 = filters.distance(h_2_pos, h_1_pos)
                        else:
                            pass
                    else:
                        pass

                    if temp_distance_1 > h_distance_1:
                        break
                    else:
                        pass

                        h_distance_1 = temp_distance_1
                        temp_distance_1 = 0.

                    counter_1 += 1

                temp_distance = 0.
                while(h_distance_2 > temp_distance_2):
                    h_2_pos = filters.rotate(h_2_pos, axis_2, shift_2, step_size_2)
                    c_2_pos = filters.rotate(c_2_pos, axis_1, shift_1, step_size_1)
                    temp_distance_2 = filters.distance(h_1_pos, h_2_pos)

                    if counter_2 == 0:
                        if temp_distance > h_distance:
                            step_size_2 = -step_size
                            h_2_pos = filters.rotate(h_2_pos, axis_2, shift_2, step_size_2)
                            c_2_pos = filters.rotate(c_2_pos, axis_1, shift_1, step_size_1)
                            temp_distance = filters.distance(h_1_pos, h_2_pos)
                        else:
                            pass
                    else:
                        pass

                    if temp_distance_2 > h_distance_2:
                        break
                    else:
                        pass

                        h_distance_2 = temp_distance_2
                        temp_distance_2 = 0.

                    counter_2 += 1

                if h_distance_1 <= h_distance_2:
                    temp_distance = h_distance_1
                else:
                    temp_distance = h_distance_2

                if temp_distance < h_distance:
                    h_distance = temp_distance
                else:
                    break

            hc_vector_1 = h_1_pos - c_1_pos
            hc_vector_2 = h_2_pos - c_2_pos
            cc_vector = c_1_pos - c_2_pos

            angle1 = filters.angle(-cc_vector, hc_vector_1)
            angle2 = filters.angle(cc_vector, hc_vector_2)

            if angle1 > angle_limit or angle2 > angle_limit:
                continue
            else:
                pass

            if h_distance >= upper_limit*length or h_distance <= lower_limit*length:
                continue

            output.append("HIT: Site {0} and Site {1} of MOF {2} allow cross-linking length of {3} and rotation steps of {4} and {5}!".format(candidates[i][0], candidates[j][0], ntpath.basename(filename), length, counter_1, counter_2))

    return output
