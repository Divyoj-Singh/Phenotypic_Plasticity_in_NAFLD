### This code genererates the topo files for the randomized networks given a list of interations.

## importing python packages: 
import itertools

## path to the output topo files: 
path_to_topo_file="./../topo_files/"


################################################################################################################
## Defining the arrays which have all the genes which are filled by reading the gene list file: 
source_array=[]
target_array=[]

## Reading the gene_list file:
with open ("./gene_list.txt") as file:
    for line in file:

        arr=line[:-1].split("\t")
        source_array.append(arr[0])
        target_array.append(arr[1])
#################################################################################################################
# Get all permutations of the interactions:
# You will have to give the interactions in any order in this,
# just the number of activation and inhibitions should be same as the core network.
array_of_interaction=[1,1,1,1,2,2,1,1,1,1]
perm = itertools.permutations([1,1,1,1,2,2,1,1,1,1],10) ## interaction array and the number of interactions:

# removing duplicates: 
set_perm=set(perm)
## now set_perm has all 45 arrays which have 1s and 2s in diffrent postions
#################################################################################################################
       
## Function to write the topo file once you give a soource array, target array and interaction array. 
def function_to_write_topo_file(src_arr,tgt_arr,intr_arr,path_to_topo_file):
   
    # Array of source and target genes for the inhibition interactions:
    inh_src=[]
    inh_tgt=[]
    
    ## finding the the inhibition interaction so that the name of the file can be named according to it. 
    for i in range(len(intr_arr)):
        if intr_arr[i]==2:
            inh_src.append(src_arr[i])
            inh_tgt.append(tgt_arr[i])
    name_of_the_file= path_to_topo_file + inh_src[0] + "_to_" + inh_tgt[0] + "_and_" + inh_src[1] + "_to_" + inh_tgt[1] + '.topo'  
    
    ## finally writing the topo file:
    with open(name_of_the_file,'w+') as file:
        file.write("Source"+'\t'+"Target"+'\t'+"Type"+'\n')
        for i in range(len(intr_arr)):
            file.write(src_arr[i]+'\t'+tgt_arr[i]+'\t'+str(intr_arr[i])+'\n')
    pass             


#################################################################################################################
### Here we are calling the above function and itterating over the interaction arrays created.
for i in set_perm:
    function_to_write_topo_file(source_array,target_array,i,path_to_topo_file)            
#################################################################################################################
