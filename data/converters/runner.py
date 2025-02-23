import os
import pathlib
import time

DIR = "/home/yoogi/GPU_Programming/community_detection/data/"

folder = pathlib.Path(DIR+"raw/")
txt_files = [str(i) for i in list(folder.rglob("*.txt"))]
csv_files = [str(i) for i in list(folder.rglob("*.csv"))]

log_file = open(DIR + "converters/log.txt", "w")

os.system("g++ " + DIR + "converters/csv_to_csr.cpp -o " + DIR + "converters/csv.out")
os.system("g++ " + DIR + "converters/txt_to_csr.cpp -o " + DIR + "converters/txt.out")

for file in txt_files:
    path = file
    name = file.split('/')[-1][:-4:]
    
    t0 = time.time()
    print(DIR + "converters/txt.out " + path + " " + DIR + "csr/" + name + ".csr")
    
    os.system(DIR + "converters/txt.out " + path + " " + DIR + "csr/" + name + ".csr")

    t1 = time.time()

    log_file.write(file.split('/')[-1] + "\t\t\t\t\t" + str(t1-t0) + "\n")
    print(name)

for file in csv_files:
    path = file
    name = file.split('/')[-1][:-4:]
    
    t0 = time.time()
    print(DIR + "converters/csv.out " + path + " " + DIR + "csr/" + name + ".csr")
    
    os.system(DIR + "converters/csv.out " + path + " " + DIR + "csr/" + name + ".csr")

    t1 = time.time()

    log_file.write(file.split('/')[-1] + "\t\t\t\t\t" + str(t1-t0) + "\n")
    print(name)



    





