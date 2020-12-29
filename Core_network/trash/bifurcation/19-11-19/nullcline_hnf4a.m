clear all
%% Species Section: All the species and an array which contains the name:
% X1: HNF4A, X2: HNF1A, X3:PPARG, X4:SREBF1,
species_array=["HNF4A","HNF1A","PPARG","SREBF1"];


%% Time Domain: 
% mention the time domain for which u want to run the ODE. 
domain = [0 1000];
bifurcation_parameter=0.00175;
%% Starting the loop for different inital conditions:
HNF4A=[];
PPARG=[];

for i=0:0.0001:2
%% nullcline wrt d(hn4a)dt=0:
IHNF4A =i ;  
IHNF1A =8;
IPPARG =3;
ISREBF1=90;
%% Calling ODE function:
[t, x] = ode45(@(t,x) hnf4a_interactions(t,x,bifurcation_parameter),domain,[IHNF4A;IHNF1A;IPPARG;ISREBF1]);

%% nullcline wrt d(hn4a)dt=0:
end_point_srebp=x(end,4);
end_point=x(end,3);
PPARG=[PPARG,end_point];
HNF4A=[HNF4A,i];
end
% PPARG=log2(PPARG);
% HNF4A=log2(HNF4A);
%% Plotting:
plot(HNF4A,PPARG) 
%% plotting
figure(1)
xlabel('HNF4A')
ylabel('PPARG')
hold on