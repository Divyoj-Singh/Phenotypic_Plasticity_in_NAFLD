clear all
%% Species Section: All the species and an array which contains the name:
% X1: HNF4A, X2: HNF1A, X3:PPARG, X4:SREBF1,
species_array=["HNF4A","HNF1A","PPARG","SREBF1"];
% The species which u want to plot: 
species_to_be_plotted=3;

%% reading the parameters file:
% change the parameter file path here:
s= tdfread('./parameters.txt','\t');

%% Parameter Set Specification:
index_no=3  %% Specify the parameter uid for which u want to plot the dynamics.
parameter_uid=s.uid(index_no)
no_of_states=s.num_states(index_no); %% as predicted by racipe.

%% Time Domain: 
% mention the time domain for which u want to run the ODE. 
domain = [0 100];
high=0;
low=0;

%% Starting the loop for different inital conditions:
for i=1:100
   
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

%% Checking if the steady state of the species is high or low:
if y(end,species_to_be_plotted)>0
    high=high+1;
elseif y(end,species_to_be_plotted)<0
    low=low+1;
end

%% Plotting:    
plot(t,y(:,species_to_be_plotted))%,'linewidth',2) 
hold on
end
figure(1)
fig_name= no_of_states;
title(fig_name)
%legend(species_array(species_to_be_plotted))
%xlim([30,1000])
xlabel('time (hour)')
ylabel('PPARG concerntration (log scale)')
% saveas(gcf, 'test.png');
% export_fig test4.png -r2500


%% Printing the proportion of high and low:
fraction_of_high=high/100
fraction_of_low=low/100

