import numpy as np


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


def common_member(candidate_1, candidate_2, candidate_3 = (-1, -1, -1, -1), candidate_4 = (-2, -2, -2, -2)):
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
