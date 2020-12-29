# this code plots the mean and std after reading from the file:
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import distance

##########################################################################################
mean_array=[]
std_array=[]
filename="bar_plot_data2.txt"
positions=[]
line_no=0
with open (filename) as f1:
    for line in f1:
        line_no+=1
        a=line[:-1].split("\t")
        mean_array.append(float(a[1]))
        std_array.append(float(a[2]))
        positions.append(line_no)


plt.bar(positions,mean_array,width=0.8, yerr=std_array,alpha=0.5, ecolor='black', capsize=10)
plt.xticks([],[])
plt.savefig('./integrtion_method.png',dpi=1000)
############################################################################################   
mean_array=[]
std_array=[]
filename="bar_plot_data1.txt"
positions=[]
line_no=0
with open (filename) as f1:
    for line in f1:
        line_no+=1
        a=line[:-1].split("\t")
        mean_array.append(float(a[1]))
        std_array.append(float(a[2]))
        positions.append(line_no)


plt.bar(positions,mean_array,width=0.8, yerr=std_array,alpha=0.5, ecolor='black', capsize=10)
plt.xticks([],[])
plt.savefig('./initial_conditions.png',dpi=1000)
###########################################################################################