clear all
%% Species Section: All the species and an array which contains the name:
% X1: HNF4A, X2: HNF1A, X3:PPARG, X4:SREBF1,
species_array=["HNF4A","HNF1A","PPARG","SREBF1"];
% The species which u want to plot: 
species_to_be_plotted=4;


%% Time Domain: 
% mention the time domain for which u want to run the ODE. 
domain = [0 1000];
bifurcation_parameter=0.65;
%% Starting the loop for different inital conditions:
for i=1:1000
% picking random initial condition for the species:    
% here we picked a random number in the range of 1-20 and converted it to log2 scale.(as done in RACIPE) 
IHNF4A = 2^(10*randn(1));
IHNF1A = 2^(10*randn(1));
IPPARG = 2^(10*randn(1));
ISREBF1= 2^(10*randn(1));

%% Calling ODE function:
[t, x] = ode15s(@(t,x) interactions(t,x,bifurcation_parameter),domain,[IHNF4A;IHNF1A;IPPARG;ISREBF1]);

% now converting the values back to log 
y=log2(x);
%% Plotting:
% plot(t,x(:,1))
% hold on
% %plot(t,x(:,2)) 
% hold on

plot(t,x(:,1)) 
hold on
%plot(t,x(:,4)) 

end
 
%% plotting
figure(1)
xlim([300,1000])
legend("PPARG")
xlabel('time')
ylabel('species concerntration')