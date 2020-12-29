clear all
%% Species Section: All the species and an array which contains the name:
% X1: HNF4A, X2: HNF1A, X3:PPARG, X4:SREBF1,
species_array=["HNF4A","HNF1A","PPARG","SREBF1"];
% The species which u want to plot: 
species_to_be_plotted=3;% 3 will correspond to PPARG


%% reading the parameters file:
% change the parameter file path here:
s= tdfread('./../parameter_files/run_1/state_3_HH_HL_LH_parameter_subset.txt','\t');

%% Parameter Set Specification:
index_no=8 %% Specify the parameter uid for which u want to plot the dynamics.
parameter_uid=s.uid(index_no)
no_of_states=s.num_states(index_no); %% as predicted by racipe.

%% Time Domain: 
% mention the time domain for which u want to run the ODE. 
domain = [0 100];

%% Starting the loop for different inital conditions:
for i=1:3
% picking random initial condition for the species:    
% here we picked a random number in the range of 1-20 and converted it to log2 scale.(as done in RACIPE) 
IHNF4A = 2^(20*rand(1));
IHNF1A = 2^(20*rand(1));
IPPARG = 2^(20*rand(1));
ISREBF1= 2^(20*rand(1));

%% Calling ODE function:
[t, x] = ode45(@(t,x) interactions(t,x,s,index_no),domain,[IHNF4A;IHNF1A;IPPARG;ISREBF1]);

% now converting the values back to log 
y=log2(x);
plot(t,y(:,species_to_be_plotted),'linewidth',2) 
hold on
end

%% Plotting and saving the figure:
figure(1)
fig_name= no_of_states;
title(fig_name)
legend(species_array(species_to_be_plotted))
%xlim([300,1000])
xlabel('time (hour)')
ylabel('PPARg concerntration (log scale)')
saveas(gcf, 'test.png');
export_fig test4.png -r2500
