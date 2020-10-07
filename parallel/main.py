import os
import numpy as np
import time
import datetime
import multiprocessing
import math
from pathlib import Path
import screening


database_path = "path"
length = 1.0
arm_length = 1.0  # only if cross-linker is symmetrical
lengths = np.array([arm_length, arm_length, arm_length, arm_length])
upper_limit = 1.1
lower_limit = 0.9
step_size = 100
angle_limit = 20.0
neighbor_distance = 1.55
n_mofs = len([name for name in os.listdir(database_path) if name.endswith(".cif")])
mof = 0
startTime = time.time()
processes = 10
my_output = "multibond_screening_" + datetime.datetime.now().strftime("%I:%M:%S-%d%m%Y") + ".log"

output_path = Path(cwd = os.getcwd())
output_path = output_path.parent



def mp_worker(queue, file_names):
    output = []
    for file_name in file_names:
        result = screening.screening(file_name, neighbor_distance, lengths, length, arm_length, upper_limit, lower_limit, step_size, angle_limit, mof, n_mofs)
        output.append(result)
    queue.put(output)


def mp_screening(file_names, processes):
    queue = multiprocessing.Queue()
    chunks = int(math.ceil(len(file_names)) / processes)
    procs = []
    for i in range(processes):
        proc = multiprocessing.Process(target=mp_worker, args=(queue, file_names[chunks*i:chunks*(i+1)]))
        procs.append(proc)
        proc.start()

    output = []
    for i in range(processes):
        output.append(queue.get())

    for i in procs:
        i.join()

    return output


def main():
    cif_files = []
    for filename in os.listdir(database_path):
        if filename.endswith(".cif"):
            file_path = os.path.join(database_path, filename)
            cif_files.append(file_path)
        else:
            continue

    temp_output = mp_screening(cif_files, processes)

    output = []
    for a in temp_output:
        for b in a:
            for c in b:
                output.append(c)


    with open(os.path.join(output_path, my_output), 'w') as f:
        for out in output:
            f.write("%s\n" % out)


if __name__ == "__main__":
        main()

print("This screening took {0} seconds to execute!".format(time.time() - startTime))
