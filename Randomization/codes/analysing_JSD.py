import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import distance

d = {}
ref_file = "./../RACIPE_output/HNF1A_to_PPARG_and_SREBF1_to_HNF4A/"
ref_distribution = []
with open(ref_file+"run_1/state_freq/JSD_distribution.txt") as f:
    for line in f:
        a = line[:-1].split("\t")
        ref_distribution.append(float(a[1]))
        
ref_plasticity_score = 0
with open(ref_file+"run_1/state_freq/states_1.txt") as f:
    for line in f:
        a = line[:-1].split("\t")
        ref_plasticity_score = 1 - (float(a[-1])/10000)
        break
print(ref_plasticity_score,(float(a[-1])/10000))
        
plasticity_score = []
for folder in os.listdir("./../RACIPE_output/"):
    if folder == "HNF1A_to_PPARG_and_SREBF1_to_HNF4A":
        continue
    d[folder] = []
    with open("./../RACIPE_output/"+folder+"/run_1/state_freq/JSD_distribution.txt") as f:
        for line in f:
            a = line[:-1].split("\t")
            d[folder].append(float(a[1]))
    with open("./../RACIPE_output/"+folder+"/run_1/state_freq/states_1.txt") as f:
        for line in f:
            a = line[:-1].split("\t")
            plasticity_score.append(1 - (float(a[-1])/10000))
            break


c = 0
jsd_dist = []
for f in d:
    c += 1
    jsd_dist.append(distance.jensenshannon(ref_distribution, d[f]))
    
ax = sns.kdeplot(jsd_dist)
plt.title("JSD distribution from the core reference network",fontname="Arial", fontsize=34)
plt.xlabel("JSD from the core network",fontname="Arial", fontsize=32)
plt.ylabel("Kernel Density Estimate",fontname="Arial", fontsize=32)
plt.xticks(fontname="Arial",fontsize=30, rotation=0)
plt.yticks(fontname="Arial",fontsize=30, rotation=0)
plt.show()
ax = sns.kdeplot(plasticity_score)
plt.title("Plasticity Score Distribution of the Randomised Networks",fontname="Arial", fontsize=34)
plt.xlabel("Plasticity Score",fontname="Arial", fontsize=32)
plt.ylabel("Kernel Density Estimate",fontname="Arial", fontsize=32)
plt.xticks(fontname="Arial",fontsize=30, rotation=0)
plt.yticks(fontname="Arial",fontsize=30, rotation=0)
plt.show()

print("jsd_dist: ",max(jsd_dist),min(jsd_dist))
print("Plasticity score: ",max(plasticity_score),min(plasticity_score))