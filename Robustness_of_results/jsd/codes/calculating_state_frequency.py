## this code calculates the frequency of HH,HL,LH and LL for each condition and print its mean and standard deviation over 3 runs.
import os
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
from textwrap import wrap
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering

filename = "./../RACIPE_output/HNF1A_to_PPARG_and_SREBF1_to_HNF4A/run_1/state_freq/states_1.txt"
path_to_RACIPE_output= "./../RACIPE_output/"
dic_states={"HH":0,"HL":0,"LH":0,"LL":0}


for network_name in os.listdir(path_to_RACIPE_output):

    dic_states={"HH":[0,0,0],"HL":[0,0,0],"LH":[0,0,0],"LL":[0,0,0]}

    for i in range(1,4):

    
        filename=path_to_RACIPE_output + network_name + "/run_"+str(i)+"/state_freq/states_1.txt"
        with open (filename) as f1:
            for line in f1:
                a=line[:-1].split("\t")
                count =int(a[-2])/10000
                dic_states[a[0]][i-1]+=count

        filename=path_to_RACIPE_output + network_name + "/run_"+str(i)+"/state_freq/states_2.txt"
        with open (filename) as f1:
            for line in f1:
                a=line[:-1].split("\t")
                count =int(a[-2])/10000
                dic_states[a[0]][i-1]+=count
                dic_states[a[1]][i-1]+=count


        filename=path_to_RACIPE_output + network_name + "/run_"+str(i)+"/state_freq/states_3.txt"
        with open (filename) as f1:
            for line in f1:
                a=line[:-1].split("\t")
                count =int(a[-2])/10000
                dic_states[a[0]][i-1]+=count
                dic_states[a[1]][i-1]+=count
                dic_states[a[2]][i-1]+=count

        filename=path_to_RACIPE_output + network_name + "/run_"+str(i)+"/state_freq/states_4.txt"
        with open (filename) as f1:
            for line in f1:
                a=line[:-1].split("\t")
                count =int(a[-2])/10000
                dic_states[a[0]][i-1]+=count
                dic_states[a[1]][i-1]+=count
                dic_states[a[2]][i-1]+=count
                dic_states[a[3]][i-1]+=count

    print(network_name)
    for i in dic_states:
        print(i,np.mean(dic_states[i]),np.std(dic_states[i]))
            