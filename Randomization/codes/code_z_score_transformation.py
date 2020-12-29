import os
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
from textwrap import wrap
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema
from sklearn.cluster import AgglomerativeClustering

# function that just normalises a given dictionary
def normalize(d, target=1.0):
    raw = sum(d.values())
    factor = target/raw
    return {key:value*factor for key,value in d.items()}

# function that returns the order of nodes and the list of parameters
def return_the_gene_and_parameter_list(path_to_prs_file):

    genes = []
    parameters = []
        
    with open(path_to_prs_file) as file:
        next(file)
        for line in file:
            a = line.split("\t")
            if "Prod_of_" in a[0]:
                genes.append(a[0][8:])
            parameters.append(a[0])

    return genes, parameters
    pass

# function to read all the data for all the state solutions combined and then plot the expression density for each gene, scatter plots for HNF4A and PPARG with or without clustering results
def collating_all_runs_for_z_score_calculation(path_to_dat_files,network_name,num_clusters,genes):

    # initializing the state dataframe that stores the expression data for all(1 to 4 state solutions) irrespective of them being in a mono, bi, tri or tetra stable cases
    state_dataframe = pd.DataFrame(columns = genes)

    # iterating over all the files with mono, bi, tri and tetra stable cases
    for i in range(1,5):
        
        # reading each file separately, getting the sub data structures and appending it to the state dataframe
        data = pd.read_csv(path_to_dat_files+network_name+"_solution_"+str(i)+".dat",delimiter="\t",header=None)
        
        for j in range(0,i):
            sub_dataframe = data[data.columns[len(genes)*j+2:len(genes)*j+(len(genes)+2)]]
            sub_dataframe.columns = genes
            state_dataframe = state_dataframe.append(sub_dataframe,ignore_index = True)
    
    mean_of_the_columns = state_dataframe.mean(axis = 0)
    std_of_the_columns = state_dataframe.std(axis = 0)

    return mean_of_the_columns, std_of_the_columns
    pass
        
# function that creates a z_normalised dat file
def creating_z_normalised_dat_files(path_to_dat_files,network_name,num_clusters,genes,mean_of_the_columns,std_of_the_columns,path_to_output_z_norm):

    # iterating over all the files with mono, bi, tri and tetra stable cases
    for i in range(1,5):
        
        # reading each file separately, getting the sub data structures and appending it to the state dataframe
        data = pd.read_csv(path_to_dat_files+network_name+"_solution_"+str(i)+".dat",delimiter="\t",header=None)
        
        z_normalised_dataframe = pd.DataFrame(columns = sorted(genes))
        z_normalised_dataframe = data[data.columns[0:2]]
        z_normalised_dataframe.columns = ["id","states"]
        
        for j in range(0,i):
            sub_dataframe = data[data.columns[len(genes)*j+2:len(genes)*j+(len(genes)+2)]]
            sub_dataframe.columns = genes
            for gene_name in genes:
                sub_dataframe[gene_name] = (sub_dataframe[gene_name] - mean_of_the_columns[gene_name])/std_of_the_columns[gene_name]
            sub_dataframe = sub_dataframe[sorted(genes)]
            
            z_normalised_dataframe = pd.concat([z_normalised_dataframe, sub_dataframe],axis=1)

        z_normalised_dataframe.to_csv(path_to_output_z_norm+network_name+"_solution_"+str(i)+"_z_norm.dat",sep="\t", header=False, index=False) 
    pass

def find_minima_in_distribution(gene_values,gene,overall_range,minima_range,accuracy):
    
    density = gaussian_kde(gene_values)
    
    data = []
    
    for i in range(overall_range[0]*accuracy,(overall_range[1]*accuracy)+1):
        data.append(density(i/accuracy)[0])
        
    lb_index = int(abs((overall_range[0]*accuracy)-(minima_range[0]*accuracy)))
    len_minima_range = int((minima_range[1]-minima_range[0])*accuracy)
    
    minima_index = argrelextrema(np.array(data[lb_index:lb_index+len_minima_range]), np.less)
    
    threshold_value = ((minima_index[0]/accuracy)+minima_range[0])[0]
    #print(overall_range,minima_range,overall_range[0]*accuracy,(overall_range[1]*accuracy)+1,lb_index,lb_index+len_minima_range,data[lb_index],data[lb_index+len_minima_range])
    
    return gene,threshold_value
    
    pass

# function to plot and store the scatter plots and the histograms for the z-normalised dat files
def plot_scatter_histograms(path_to_output_z_norm,network_name,num_clusters,genes,run,path_to_plots,accuracy,thresholds,overall_range,minima_range):

    # initializing the state dataframe that stores the expression data for all(1 to 4 state solutions) irrespective of them being in a mono, bi, tri or tetra stable cases
    state_dataframe = pd.DataFrame(columns = sorted(genes))

    # iterating over all the files with mono, bi, tri and tetra stable cases
    for i in range(1,5):
        # reading each file separately, getting the sub data structures and appending it to the state dataframe
        data = pd.read_csv(path_to_output_z_norm+network_name+"_solution_"+str(i)+"_z_norm.dat",delimiter="\t",header=None)
        for j in range(0,i):
            sub_dataframe = data[data.columns[len(genes)*j+2:len(genes)*j+(len(genes)+2)]]
            sub_dataframe.columns = sorted(genes)
            state_dataframe = state_dataframe.append(sub_dataframe, ignore_index = True)

    #print("Now running on: ",network_name, ";  Iteration number: ",run+1)
    
    s = network_name+"\t"
    t = network_name+"\t"
    for g in genes:
        #gene,threshold_value = #find_minima_in_distribution(state_dataframe[g],g,overall_range,minima_range,accuracy)
        #s += str(threshold_value)+"\t"
        #t += str(gene)+"\t"
        #print(threshold_value)
        thresholds[network_name][g] = 0#threshold_value
    #print(t)
    #print(s)

    # plotting the distribution plots for each of the genes
    ax = sns.distplot(state_dataframe["HNF4A"],hist=False,kde_kws={"color": "r", "label": "HNF4A"})
    ax = sns.distplot(state_dataframe["PPARG"],hist=False,kde_kws={"color": "b", "label": "PPARG"})
    ax = sns.distplot(state_dataframe["HNF1A"],hist=False,kde_kws={"color": "y", "label": "HNF1A"})
    ax = sns.distplot(state_dataframe["SREBF1"],hist=False,kde_kws={"color": "g", "label": "SREBF1"})
    plt.title(network_name+"__Replicate_"+str(run+1))
    plt.xlabel("Z-normalised Expression")
    plt.ylabel("Frequency")
    plt.savefig(path_to_plots+network_name+"__Replicate_"+str(run+1)+"_histogram.png",dpi=1000)
    plt.close()
            
    """for pair in [["HNF4A","PPARG"],["HNF1A","PPARG"],["HNF4A","SREBF1"],["HNF1A","SREBF1"],["PPARG","SREBF1"],["HNF4A","HNF1A"]]:
        sc = plt.scatter(state_dataframe[pair[0]],state_dataframe[pair[1]],marker="o",s=1)
        plt.title(network_name+"__Replicate_"+str(run+1))
        plt.xlabel(pair[0])
        plt.ylabel(pair[1])
        plt.ylim(-3, 2)
        plt.xlim(-3, 2)
        plt.savefig(path_to_plots+network_name+"_Replicate_"+str(run+1)+"_"+pair[0]+"_"+pair[1]+"_scatter.png",dpi=1000)
        plt.close()"""

    # plotting the scatter plot with clustering level data on the HNF4A and PPARG space
    """cluster = AgglomerativeClustering(n_clusters=num_clusters, affinity='euclidean', linkage='ward')
    cluster.fit_predict(state_dataframe)
    plt.scatter(state_dataframe["HNF4A"],state_dataframe["PPARG"], c=cluster.labels_, cmap='rainbow',marker="o",s=0.5)
    plt.title(network_name+"__Replicate_"+str(run+1))
    plt.xlabel("HNF4A")
    plt.ylabel("PPARG")
    plt.savefig(path_to_plots+network_name+"__Replicate_"+str(run+1)+"_clusters.png")
    plt.close()"""
    
    return thresholds
    
    pass

# generating the frequencies for each combination of states
def generating_state_frequencies(path_to_output_z_norm,network_name,genes,thresholds,states,genes_to_probe,path_to_state_frequencies):

    
    
    list_of_indexs = sorted([genes.index(genes_to_probe[0]),genes.index(genes_to_probe[1])])
    
    threshold_array_2_gene = [thresholds[network_name][genes[list_of_indexs[0]]],thresholds[network_name][genes[list_of_indexs[1]]]]
    
    state_dataframe = pd.DataFrame(columns = sorted(genes))
    
    for i in range(1,5):
    
        combinations  = itertools.combinations_with_replacement(states,i)
        
        dict_freq = {}

        for j in combinations:
            dict_freq[j] = 0  
            
        data = pd.read_csv(path_to_output_z_norm+network_name+"_solution_"+str(i)+"_z_norm.dat",delimiter="\t",header=None)
        index_array = sorted(list(range(2+list_of_indexs[0],i*len(genes)+2,len(genes))) + list(range(2+list_of_indexs[1],i*len(genes)+2,len(genes))))
        threshold_array = np.tile(threshold_array_2_gene,i)
       
        for index, row in data.iterrows():
        
            expression_array = []
            
            for j in index_array:
                expression_array.append(row[j])
            
            condition_array = np.array(expression_array) >= threshold_array
            
            condition_array = condition_array.astype('str')
            
            condition_array[condition_array == "True"] = "H"
            
            condition_array[condition_array == "False"] = "L"
            
            condition_string = ""
            condition_string = condition_string.join(condition_array)
            condition = wrap(condition_string,2)
            
            for j in dict_freq:
                l = list(j)
                l.sort()
                condition.sort()
                if l == condition:
                    dict_freq[j] += 1
                    pass
            
        total_freq = sum(dict_freq.values())
        dict_freq_norm = normalize(dict_freq)
        with open (path_to_state_frequencies+"/states_"+str(i)+".txt","w") as f:
            for j in dict_freq:
                for k in j:
                    f.write(k+"\t")
                f.write(str(round(dict_freq_norm[j],5))+"\t"+str(dict_freq[j])+"\t"+str(total_freq)+"\n")
    pass

# this function reads the parameters file 
def reading_the_parameters_file(path_to_prs_file):
    parameter_data = pd.read_csv(path_to_prs_file,delimiter="\t",header=None)
    return parameter_data

# generating the frequencies for each combination of states
def generating_the_parameters_file(path_to_output_z_norm,network_name,thresholds,genes,states,parameter_data,path_to_relative_stab_params):
    
    list_of_indexs = sorted([genes.index(genes_to_probe[0]),genes.index(genes_to_probe[1])])
    
    threshold_array_2_gene = [thresholds[network_name][genes[list_of_indexs[0]]],thresholds[network_name][genes[list_of_indexs[1]]]]
    
    state_dataframe = pd.DataFrame(columns = sorted(genes))
    
    for i in range(1,5):
    
        combinations  = itertools.combinations_with_replacement(states,i)
    
        dict_freq = {}

        for j in combinations:
            dict_freq[j] = 0    
            
        data = pd.read_csv(path_to_output_z_norm+network_name+"_solution_"+str(i)+"_z_norm.dat",delimiter="\t",header=None)
        index_array = sorted(list(range(2+list_of_indexs[0],i*len(genes)+2,len(genes))) + list(range(2+list_of_indexs[1],i*len(genes)+2,len(genes))))
        threshold_array = np.tile(threshold_array_2_gene,i)
       
        for index, row in data.iterrows():
        
            expression_array = []
            
            for j in index_array:
                expression_array.append(row[j])
            
            condition_array = np.array(expression_array) >= threshold_array
            
            condition_array = condition_array.astype('str')
            
            condition_array[condition_array == "True"] = "H"
            
            condition_array[condition_array == "False"] = "L"
            
            condition_string = ""
            condition_string = condition_string.join(condition_array)
            condition = wrap(condition_string,2)
            
            for j in dict_freq:
                l = list(j)
                l.sort()
                condition.sort()
                if l == condition:
                    identifier = ""
                    for k in j:
                        identifier += k + "_"
                    parameter_uid = int(row[0])
                    with open(path_to_relative_stab_params+"state_"+str(i)+"_"+identifier+"parameter_subset.txt","a+") as f:
                        string_to_be_written = ""
                        for j in parameter_data.loc[parameter_uid-1]:
                            string_to_be_written += str(j) + "\t"
                        f.write(string_to_be_written[:-1]+"\n")


path_to_RACIPE_output = "./../RACIPE_output/"
num_clusters = 4
number_of_runs = 3
accuracy = 1000
states = ["HH","HL","LH","LL"]
overall_range = [-3,3]
minima_range = [0,1]
thresholds = {}
genes_to_probe = ["HNF4A","PPARG"]
list_of_files_to_store = ["state_2_HH_HL_parameter_subset","state_2_HH_LH_parameter_subset","state_2_HL_LH_parameter_subset","state_2_HL_LL_parameter_subset","state_2_LH_LL_parameter_subset","state_3_HH_HL_LH_parameter_subset","state_3_HL_LH_LL_parameter_subset","state_4_HH_HL_LH_LL_parameter_subset"]

for network_name in os.listdir(path_to_RACIPE_output):

    #print("Running Z-score transformation on the network: ",network_name)
    
    thresholds[network_name] = {}
    
    for run in range(number_of_runs):
        
        path_to_dat_files = path_to_RACIPE_output+network_name+"/"+"run_"+str(run+1)+"/"
        path_to_output_z_norm = path_to_RACIPE_output+network_name+"/"+"run_"+str(run+1)+"/"+"Z_norm/"
        path_to_plots = path_to_RACIPE_output+network_name+"/"+"run_"+str(run+1)+"/"+"plots/"
        path_to_state_frequencies = path_to_RACIPE_output+network_name+"/"+"run_"+str(run+1)+"/"+"state_freq/"
        path_to_relative_stab_params = path_to_RACIPE_output+network_name+"/"+"run_"+str(run+1)+"/"+"rel_stab_params/"
        path_to_relative_stab_params_selc = path_to_RACIPE_output+network_name+"/"+"run_"+str(run+1)+"/"+"rel_stab_params_selc/"
        path_to_prs_file = path_to_RACIPE_output+network_name+"/"+"run_"+str(run+1)+"/"+network_name+".prs"
        path_to_parameters_file = path_to_RACIPE_output+network_name+"/"+"run_"+str(run+1)+"/"+network_name+"_parameters.dat"
        
        if not os.path.exists(path_to_output_z_norm):
            os.makedirs(path_to_output_z_norm)
        if not os.path.exists(path_to_plots):
            os.makedirs(path_to_plots)
        if not os.path.exists(path_to_state_frequencies):
            os.makedirs(path_to_state_frequencies)
        if not os.path.exists(path_to_relative_stab_params):
            os.makedirs(path_to_relative_stab_params)
        if not os.path.exists(path_to_relative_stab_params_selc):
            os.makedirs(path_to_relative_stab_params_selc)
            
        genes,parameters = return_the_gene_and_parameter_list(path_to_prs_file)
                
        header_row = "uid" + "\t" +"num_states" + "\t"
        for i in parameters:
            header_row += i + "\t"
        header_row = header_row[:-1]+"\n"
        
        mean_of_the_columns,std_of_the_columns = collating_all_runs_for_z_score_calculation(path_to_dat_files,network_name,num_clusters,genes)
        
        creating_z_normalised_dat_files(path_to_dat_files,network_name,num_clusters,genes,mean_of_the_columns,std_of_the_columns,path_to_output_z_norm)
        
        thresholds = plot_scatter_histograms(path_to_output_z_norm,network_name,num_clusters,genes,run,path_to_plots,accuracy,thresholds,overall_range,minima_range)
        
        generating_state_frequencies(path_to_output_z_norm,network_name,genes,thresholds,states,genes_to_probe,path_to_state_frequencies)
        
        """parameter_data = reading_the_parameters_file(path_to_parameters_file)
        
        generating_the_parameters_file(path_to_output_z_norm,network_name,thresholds,genes,states,parameter_data,path_to_relative_stab_params)
            
        for filename in os.listdir(path_to_relative_stab_params):
            if "parameter_subset" in filename and filename[:-4] in list_of_files_to_store:
            
                with open(path_to_relative_stab_params_selc+filename,"w") as file2:
                    file2.write(header_row)
                    with open(path_to_relative_stab_params+filename) as file1:
                        for line in file1:
                            file2.write(line)"""