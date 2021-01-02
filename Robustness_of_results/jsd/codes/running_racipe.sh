## The code needs a folder on which topo files are kept and a folder in which list of topo files is kept. It makes a folder RACIPE_output and keeps the outputs ( kde plots, selected relative stability paramters, state frequency and z normalised steady states) in it. 

num_of_parameters=10000  ## no of parametres

## Paths
path_to_codes_folder="/media/csb/D/NAFLD2/Perturbation/RACIPE/codes"
path_to_list_of_topo_files="./../list_of_topo_files/new_list.txt"
path_to_topo_file_folder="./../topo_files/"
path_to_RACIPE_output_folder="/media/csb/D/NAFLD2/jsd_calculation/racipe/RACIPE_output_init_cond_1000/"
path_to_RACIPE="/home/csb/RACIPE-1.0-master/RACIPE"

while read network_identifier; do

    # iterating over all the network instances for which RACIPE needs to be run
    #-num_paras $num_of_parameters -seed $variable
	echo Running on the network: $network_identifier
	
	path_to_topo_file=$path_to_topo_file_folder$network_identifier.topo
	

	for variable in {1..3} ### the RUN index 
		do
			path_to_run_folder=$path_to_RACIPE_output_folder$network_identifier/run_$variable/
			path_to_run_topo_file=$path_to_run_folder$network_identifier.topo
			
			mkdir -p $path_to_run_folder
			cp $path_to_topo_file $path_to_run_folder
			cd $path_to_run_folder
			
			## THE FOLLOWING LINE RUNS THE RACIPE. LOOK AT IT PROPERLY.
			$path_to_RACIPE $network_identifier.topo -num_paras $num_of_parameters -seed $variable -num_stability 4 -num_ode 1000 & 
            #-maxF 10 &
			cd $path_to_codes_folder
		done

done < $path_to_list_of_topo_files


#python3 ./code_z_score_transformation.py
