clear all
%% Species Section: All the species and an array which contains the name:
% X1: HNF4A, X2: HNF1A, X3:PPARG, X4:SREBF1,
species_array=["HNF4A","HNF1A","PPARG","SREBF1"];

%% Time Domain: 
% mention the time domain for which u want to run the ODE. 
domain = [0 100000];
bifurcation_parameter=0.38;
%% Starting the loop for different inital conditions:
%for i=1:1
% picking random initial condition for the species:    
% here we picked a random number in the range of 1-20 and converted it to log2 scale.(as done in RACIPE) 

IHNF4A = 62;  %2^(20*rand(1));
IHNF1A =   15;%2^(20*rand(1));
IPPARG = 5;%2^(20*rand(1));
ISREBF1=   90;%2^(20*rand(1));

%% Calling ODE function:
[t, x] = ode45(@(t,x) interactions(t,x),domain,[IHNF4A;IHNF1A;IPPARG;ISREBF1]);

% now converting the values back to log 
%y=log2(x);
% for j=1:(size(x,1))
%     
% % finding the ratio of PPARG(1) to HNF4A(3):      
% %x(j,5)=t(j,1);
% end

%% Plotting:
plot(t,x(:,1),'b','linewidth',2)
ylabel('HNF4A')
yyaxis right

plot(t,x(:,3),'r','linewidth',2)
ylabel('PPARG')
%% plotting
figure(1)
% t1=linspace(0,416.66);
% xticks(t1)
xlabel('time(hours)')
%saveas(gcf, 'test.png');
%export_fig test4.png -r2500