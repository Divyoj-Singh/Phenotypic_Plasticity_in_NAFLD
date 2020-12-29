clear all

th_2=0;
th_1=0;

no_of_initial_conditions=1000;

%% output file:
path_to_fractions_file='./matlab_output/fraction_for_HL_LH_HH_LL.txt';
header=["bifurcation_parameter","HH","HL","LH","LL"];
fid = fopen(path_to_fractions_file,'wt');
fprintf(fid,'%s\t',header);
fprintf(fid,'\n');
fclose(fid);


%% Running loop over different parameter:
for bifurcation_parameter=0:0.01:5
    
%% Parameter Set Specification:
index_no=bifurcation_parameter

%% Time Domain:
% mention the time domain for which u want to run the ODE. 
domain = [0 1000];
state = zeros(1,4);

%% Starting the loop for different inital conditions:
parfor j=1:no_of_initial_conditions

    % picking random initial condition for the species:    
    % here we picked a random number in the ra4nge of 1-20 and converted it to log2 scale.(as done in RACIPE) 
    IA = 2^(20*rand(1));
    IB = 2^(20*rand(1));

%% Calling ODE function:
[t, x] = ode45(@(t,x) interactions(t,x,bifurcation_parameter),domain,[IA;IB]);

% now converting the values back to log 
x=log2(x);
ind = 0;
s1 = zeros(1,4);
if (x(end,1)> th_1) && (x(end,2)>th_2)
    ind = 1;

elseif (x(end,1)> th_1) && (x(end,2)<th_2) 
    ind = 2;

elseif (x(end,1)< th_1) && (x(end,2)>th_2) 
    ind = 3;

elseif (x(end,1)< th_1) && (x(end,2)<th_2) 
    ind = 4; 
end

s1(ind) = 1;
state = [state; s1];
end
count = struct('HH',sum(state(:,1)),'HL',sum(state(:,2)),'LH',sum(state(:,3)),'LL',sum(state(:,4)));
fid = fopen(path_to_fractions_file,'a+');
printing_array=[bifurcation_parameter,(count.HH/no_of_initial_conditions),(count.HL/no_of_initial_conditions),(count.LH/no_of_initial_conditions),(count.LL/no_of_initial_conditions)];
fprintf(fid,'%g\t',printing_array);
fprintf(fid,'\n');
fclose(fid);
end

