import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import distance

d = {}
ref_file = "./../RACIPE_output/HNF1A_to_PPARG_and_SREBF1_to_HNF4A/"
ref_distribution = []
with open(ref_file+"run_2/state_freq/JSD_distribution.txt") as f:
    for line in f:
        a = line[:-1].split("\t")
        ref_distribution.append(float(a[1]))

        

for folder in os.listdir("./../RACIPE_output/"):
    d[folder] = []
    with open("./../RACIPE_output/"+folder+"/run_2/state_freq/JSD_distribution.txt") as f:
        for line in f:
            a = line[:-1].split("\t")
            d[folder].append(float(a[1]))
            
    with open("./../RACIPE_output/"+folder+"/run_2/state_freq/states_1.txt") as f:
        for line in f:
            a = line[:-1].split("\t")
            
            break


c = 0
jsd_dist = []
for f in d:
    c += 1
    jsd_dist.append(distance.jensenshannon(ref_distribution, d[f]))
    print("jsd",f,distance.jensenshannon(ref_distribution, d[f]))
    

print("jsd_dist: ",max(jsd_dist),min(jsd_dist))