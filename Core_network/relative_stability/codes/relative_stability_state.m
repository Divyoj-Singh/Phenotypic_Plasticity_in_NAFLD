clear all
%% Species Section: All the species and an array which contains the name:
% X1: HNF4A, X2: HNF1A, X3:PPARG, X4:SREBF1,
species_array=["HNF4A","HNF1A","PPARG","SREBF1"];
%% The names of runs in the following array:
run_array=["run_1"];
%% Reading the threshold file and storing thresholds(decided by RACIPE analysis):
path_to_threshold_file='./threshold.txt'
thresh = tdfread(path_to_threshold_file,'\t');

th_3=thresh.PPARG(1);
th_1=thresh.HNF4A(1);

%% Giving the parameter file and RUN number:
for k=1:size(run_array,1)
run_id=run_array(k);
path = './../parameter_files/'+run_id+'/';
files = dir (path);

L = length (files);
for j=1:L
    if files(j).isdir==1
        continue
    end
    parameter_file_name=files(j).name

%% Opening the output file and writing the header in it:
path_to_fractions_file='./../matlab_output/'+run_id+'/'+'fraction_for_'+parameter_file_name;
header=["index_no","uid","no_of_states","HH","HL","LH","LL"];
fid = fopen(path_to_fractions_file,'wt');
fprintf(fid,'%s\t',header);
fprintf(fid,'\n');
fclose(fid);

%% reading the parameters file:

% change the parameter file path here:
s= tdfread('./../parameter_files/'+run_id+'/'+parameter_file_name,'\t');

%% Running loop over different parameter:
for i=1:size(s.uid,1)
%% Parameter Set Specification:
index_no=i;
if rem(i,100)==0
	disp(index_no)
end
parameter_uid=s.uid(index_no);
no_of_states=s.num_states(index_no); %% as predicted by racipe.

%% Time Domain:
% mention the time domain for which u want to run the ODE. 
domain = [0 1000];
%% variable for storing the proportions:
state = zeros(1,4);

%% Starting the loop for different inital conditions:
no_of_initial_conditions=1000;
parfor j=1:no_of_initial_conditions
% picking random initial condition for the species:    
% here we picked a random number in the ra4nge of 1-20 and converted it to log2 scale.(as done in RACIPE) 
IHNF4A = 2^(20*randn(1));
IHNF1A  = 2^(20*randn(1));
IPPARG = 2^(20*randn(1));
ISREBF1= 2^(20*randn(1));

%% Calling ODE function:
[t, x] = ode45(@(t,x) interactions(t,x,s,index_no),domain,[IHNF4A;IHNF1A;IPPARG;ISREBF1]);

% now converting the values back to log 
x=log2(x);

%% Comparing the stready state values with the thresholds and deciding the phase:
ind = 0;
s1 = zeros(1,4);
if (x(end,1)> th_1) && (x(end,3)>th_3)
    ind = 1;

elseif (x(end,1)> th_1) && (x(end,3)<th_3) 
    ind = 2;

elseif (x(end,1)< th_1) && (x(end,3)>th_3) 
    ind = 3;

elseif (x(end,1)< th_1) && (x(end,3)<th_3) 
    ind = 4; 
end

s1(ind) = 1;
state = [state; s1];
end

%% counting and writing these relative proportions:
count = struct('HH',sum(state(:,1)),'HL',sum(state(:,2)),'LH',sum(state(:,3)),'LL',sum(state(:,4)));
fid = fopen(path_to_fractions_file,'a+');
printing_array=[index_no,parameter_uid,no_of_states,(count.HH/no_of_initial_conditions),(count.HL/no_of_initial_conditions),(count.LH/no_of_initial_conditions),(count.LL/no_of_initial_conditions)];
fprintf(fid,'%g\t',printing_array);
fprintf(fid,'\n');
fclose(fid);
end
end
end
