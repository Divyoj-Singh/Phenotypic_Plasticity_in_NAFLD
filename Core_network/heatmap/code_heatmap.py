## importing the relevant python modules:
import pandas as pd
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import numpy as np

## path to dat file:
path_to_dat_files="./run3/"

network_name="network_v2_4a_1a_g_f1"
genes=["HNF4A", "HNF1A", "PPARG" ,"SREBP1"]

def making_pd_data_frame(path_to_dat_files,network_name,genes):

    # initializing the state dataframe that stores the expression data for all(1 to 4 state solutions) irrespective of them being in a mono, bi, tri or tetra stable cases
    state_dataframe = pd.DataFrame(columns = genes)

    # iterating over all the files with mono, bi, tri and tetra stable cases
    for i in range(1,5):
        
        # reading each file separately, getting the sub data structures and appending it to the state dataframe
        ## tak care while reading files

        data = pd.read_csv(path_to_dat_files+network_name+"_solution_"+str(i)+"_z_norm.dat",delimiter="\t",header=None)
        
        for j in range(0,i):
            sub_dataframe = data[data.columns[len(genes)*j+2:len(genes)*j+(len(genes)+2)]]
            ## defining the subdata frame for the second, third, fourth solutions
            sub_dataframe.columns = genes
            state_dataframe = state_dataframe.append(sub_dataframe,ignore_index = True)

    return state_dataframe
    pass

# calling the function and sorting the dataframe:
pd_dataframe = making_pd_data_frame(path_to_dat_files,network_name,genes)
pd_dataframe = pd_dataframe.sort_values(["HNF4A", "HNF1A", "PPARG" ,"SREBP1"], ascending=[1, 1,1,1])


## plotting the heat map   
    
all_4 = sns.clustermap(pd_dataframe, metric="euclidean",cmap="seismic",robust=True,figsize=(6, 7),col_cluster=True, yticklabels=False)
# putting in data into the heat map
all_4.savefig( path_to_dat_files+"Clustermap_of_4_genes.png",dpi=2000)
# saving fig

a4_g = sns.clustermap(pd_dataframe[["HNF4A","PPARG"]], metric="euclidean",cmap="seismic",robust=True,figsize=(6, 7),col_cluster=True, yticklabels=False)
# putting in data into the heat map
a4_g.savefig( path_to_dat_files+"Clustermap_of_HNF4A_and_PPARG.png",dpi=2000)
# saving fig = sns.clustermap(pd_dataframe, metric="euclidean",cmap="seismic",robust=True,figsize=(6, 7),col_cluster=True, yticklabels=False)

a1_a4 = sns.clustermap(pd_dataframe[["HNF4A","HNF1A"]], metric="euclidean",cmap="seismic",robust=True,figsize=(6, 7),col_cluster=False, yticklabels=False)
# putting in data into the heat map
a1_a4.savefig( path_to_dat_files+"Clustermap_of_HNF4A_and_HNF1A.png",dpi=2000)
# saving fig

a4_f = sns.clustermap(pd_dataframe[["HNF4A","SREBP1"]], metric="euclidean",cmap="seismic",robust=True,figsize=(6, 7),col_cluster=False, yticklabels=False)
# putting in data into the heat map
a4_f.savefig( path_to_dat_files+"Clustermap_of_HNF4A_and_SREBP1.png",dpi=2000)
# saving fig

a4_f = sns.clustermap(pd_dataframe[["HNF1A","SREBP1"]], metric="euclidean",cmap="seismic",robust=True,figsize=(6, 7),col_cluster=False, yticklabels=False)
# putting in data into the heat map
a4_f.savefig( path_to_dat_files+"Clustermap_of_HNF1A_and_SREBP1.png",dpi=2000)
# saving fig

a4_f = sns.clustermap(pd_dataframe[["HNF1A","PPARG"]], metric="euclidean",cmap="seismic",robust=True,figsize=(6, 7),col_cluster=False, yticklabels=False)
# putting in data into the heat map
a4_f.savefig( path_to_dat_files+"Clustermap_of_HNF1A_and_PPARG.png",dpi=2000)
# saving fig

a4_f = sns.clustermap(pd_dataframe[["PPARG","SREBP1"]], metric="euclidean",cmap="seismic",robust=True,figsize=(6, 7),col_cluster=False, yticklabels=False)
# putting in data into the heat map
a4_f.savefig( path_to_dat_files+"Clustermap_of_PPARG_and_SREBP1.png",dpi=2000)
# saving fig