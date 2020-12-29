This code was used for generating the figure showing the 3 states: monostable, bistable and tristable.
This code takes parameter set  file with the following header:
uid	num_states	Prod_of_HNF4A	Prod_of_HES6	Prod_of_PPARG	Prod_of_SREBF1	Deg_of_HNF4A	Deg_of_HES6	Deg_of_PPARG	Deg_of_SREBF1	Trd_of_HNF4AToHNF4A	Num_of_HNF4AToHNF4A	Act_of_HNF4AToHNF4A	Trd_of_HES6ToHNF4A	Num_of_HES6ToHNF4A	Act_of_HES6ToHNF4A	Trd_of_SREBF1ToHNF4A	Num_of_SREBF1ToHNF4A	Inh_of_SREBF1ToHNF4A	Trd_of_HNF4AToHES6	Num_of_HNF4AToHES6	Act_of_HNF4AToHES6	Trd_of_HES6ToPPARG	Num_of_HES6ToPPARG	Inh_of_HES6ToPPARG	Trd_of_SREBF1ToPPARG	Num_of_SREBF1ToPPARG	Act_of_SREBF1ToPPARG	Trd_of_PPARGToPPARG	Num_of_PPARGToPPARG	Act_of_PPARGToPPARG	Trd_of_PPARGToSREBF1	Num_of_PPARGToSREBF1	Act_of_PPARGToSREBF1

The header will be read as variable in the matlab code.

The hill function  computes the hill function for that value for given paramters and is called by the interaction function .

The interaction fuction contains all the interactions in the network, and there parameters.

Dynamical behaviour script
The parameter set file, parameter set and species are to be specified here.
other parameters like time domain, no of initial conditions are also specified here.
The dynamical behaviour script calls the interactions fun and plots the graphs.
