import os
import pathlib

DIR = "/home/yoogi/GPU_Programming/community_detection/data/"

folder = pathlib.Path(DIR+"raw/")
txt_files = [str(i) for i in list(folder.rglob("*.txt"))]
csv_files = [str(i) for i in list(folder.rglob("*.csv"))]

os.system("g++ " + DIR + "converters/csv_to_csr.cpp -o " + DIR + "converters/csv.out")

for file in txt_files:
    path = file
    name = file.split('/')[-1][:-4:]
    





