## importing the python modules:
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import pandas as pd
import pylab as p
import os
import seaborn as sns

# function to plot the kde plot:
def plotting_the_kde_plot_for_the_given_cases(input_directory,output_directory,network_name):
    
    for filename in os.listdir(input_directory):

        print(filename)
        
        path_to_file = input_directory + filename

        with open (path_to_file) as file:
            next(file)
            HH=[]
            HL=[]
            LH=[]
            LL=[]
            for line in file:
                a=line[:-1].split("\t")
                HH.append(float(a[3]))
                HL.append(float(a[4]))
                LH.append(float(a[5]))
                LL.append(float(a[6]))

            total_freq = len(HH)
                    
            if 'HH' in filename:
                if sum(HH)!=0:
                    ax=sns.kdeplot(HH,label="HH")
            if 'HL' in filename:
                if sum(HL)!=0:
                    ax=sns.kdeplot(HL,label="HL")
            if 'LH' in filename:
                if sum(LH)!=0:
                    ax=sns.kdeplot(LH,label="LH")
            if 'LL' in filename:
                if sum(LL)!=0:
                    ax=sns.kdeplot(LL,label="LL")
            ax.set_ylim(0,3)
            plt.xlabel('proportion of steady state out of 1000 diff initial condition')
            plt.ylabel('frequency')
            ax.set_ylim(0,3)
            plt.title(filename[21:-4]+' for '+network_name+'\n'+str(total_freq))  
            plt.legend
            plt.savefig(output_directory+filename[:-4]+'.png') 
            plt.close()
    pass

network_name="network_v2_4a_1a_g_f1"
path_to_main_folder = "./../"
input_directory = path_to_main_folder+folder+"/matlab_output/run_1/"
output_directory = path_to_main_folder+folder+"/graphs_results/run_1_normal/"
plotting_the_kde_plot_for_the_given_cases(input_directory,output_directory,folder)