clear all
%% Species Section: All the species and an array which contains the name:
% X1: HNF4A, X2: HNF1A, X3:PPARG, X4:SREBF1,
species_array=["HNF4A","HNF1A","PPARG","SREBF1"];

%% reading the parameters file:
% change the parameter file path here:
s= tdfread('./../parameter_files/run_1/bistable_parameter_set_file.txt','\t');

%% Parameter Set Specification:
index_no=2%% Specify the parameter uid for which u want to plot the dynamics.
parameter_uid=s.uid(index_no)
no_of_states=s.num_states(index_no); %% as predicted by racipe.
%% Time Domain: 
% mention the time domain for which u want to run the ODE. 
domain = [0 1000];

%% Starting the loop for different inital conditions:
for i=1:100
% picking random initial condition for the species:    
% here we picked a random number in the range of 1-20 and converted it to log2 scale.(as done in RACIPE) 
IHNF4A = 2^(20*rand(1));
IHNF1A  = 2^(20*rand(1));
IPPARG = 2^(20*rand(1));
ISREBF1= 2^(20*rand(1));

%% Calling ODE function:
[t, x] = ode45(@(t,x) interactions(t,x,s,index_no),domain,[IHNF4A;IHNF1A;IPPARG;ISREBF1]);

% now converting the values back to log 
x=log2(x);

%%  Uncomment this section for plotting the ratio of PPARG to HNF4A for one parameter set for all initial condition.
for j=1:(size(x,1))
    
% finding the ratio of PPARG(1) to HNF4A(3):      
x(j,5)=x(j,3)/x(j,1);
end
% skipping the initial time since the values change a lot in that part and
% steady values are not clear.
plot(t(100:end),x(100:end,5))
hold on
% 100:end
end

%% Plotting
figure(1)
fig_name= no_of_states;
title(fig_name)
legend("Ratio of PPARG to HNF4A")
xlabel('time')
ylabel('species concerntration')