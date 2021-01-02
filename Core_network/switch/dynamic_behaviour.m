clear all
%% Species Section: All the species and an array which contains the name:
% X1: HNF4A, X2: HNF1A, X3:PPARG, X4:SREBF1,
species_array=["HNF4A","HNF1A","PPARG","SREBF1"];

%% Time Domain and other parameters: 
% mention the time domain for which u want to run the ODE. 
domain = [0 10000];

% Steady-State values of PPARG in absence of noise: 
high_state=15.0294;
low_state=1.5820;

%% Initial Condition:
IHNF4A = 62; 
IHNF1A = 15;
IPPARG = 8;
ISREBF1= 90;

%% Calling ODE function:
[t, x] = ode23(@(t,x) interactions(t,x),domain,[IHNF4A;IHNF1A;IPPARG;ISREBF1]);


%% Plotting:
plot(t,x(:,1),'b','linewidth',2)
ylabel('HNF4A')
yyaxis right
plot(t,x(:,3),'r','linewidth',2)
ylabel('PPARG')

%% Saving the plot:
saveas(gcf,filename)
close
figure(2)
t1=linspace(0,416.66);
xticks(t1)
label('time(hours)')
saveas(gcf, 'test.png');
export_fig test4.png -r2500

