import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import pandas as pd
import pylab as p
import os
import seaborn as sns

def plotting_the_kde_plot_for_the_given_cases(input_directory,output_directory):
    
    for filename in os.listdir(input_directory):

        print(filename)
        
        path_to_file = input_directory + filename

        with open (path_to_file) as file:
            next(file)
            HH=[]
            HL=[]
            LH=[]
            LL=[]
            k_pparg=[]
            for line in file:
                a=line[:-1].split("\t")
                k_pparg.append(float(a[0]))
                HH.append(float(a[1]))
                HL.append(float(a[2]))
                LH.append(float(a[3]))
                LL.append(float(a[4]))

            total_freq = len(HH)
            plt.plot(k_pparg,HH,label="HH")
            plt.plot(k_pparg,HL,label="HL")            
            plt.plot(k_pparg,LH,label="LH")
            plt.plot(k_pparg,LL,label="LL")
            plt.legend(title='states')
            plt.xlabel('k_pparg')
            plt.ylabel('proportion of steady state out of 1000 diff initial condition')
            #plt.title(filename[21:-4]+' for '+'\n'+str(total_freq))  
            plt.savefig(output_directory+filename[:-4]+'.png') 
            plt.close()
    pass

#if not os.path.exists(output_directory):
#        os.makedirs(output_directory)


input_directory ="./matlab_output/"
output_directory ="./graphs_results/"
if not os.path.exists(output_directory):
    os.makedirs(output_directory)
plotting_the_kde_plot_for_the_given_cases(input_directory,output_directory)