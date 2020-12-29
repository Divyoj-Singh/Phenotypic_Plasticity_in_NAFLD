function dxdt = interactions(t,x,bifurcation_parameter)
%% X1: HNF4A, X2: HNF1A, X3:PPARG, X4:SREBF1,
dxdt = zeros(4,1);
%%% production rate
g_HNF4A = 0.08540700422860051;
g_HNF1A = 0.10705569053323611;
g_PPARG = 0.06283939658384676;
g_SREBF1= 1.1058054056141339;

% Degradation Rate
k_HNF4A = 0.049119799999999997;   %% #% 0.0479
k_HNF1A = 0.049445600000000006; %% #% 0.05
k_PPARG = 0.0001;%0.07670399999999999;  %% #% 0.06
k_SREBF1 = 0.15984530000000001; %% #% 0.1732

%Thresholds
x_HNF4A_HNF4A = 3;%7.445485000000001;
x_SREBF1_SREBF1 = 4.6490279999999995;
x_PPARG_PPARG = 6.32011;
x_HNF1A_HNF1A = 15.796237;

x_SREBF1_HNF4A = 7.067421;
x_HNF1A_PPARG = 26.740759999999998;

x_HNF4A_HNF1A = 5.108258999999999;
x_HNF1A_HNF4A = 14.169678;
x_SREBF1_PPARG = 26.013288;
x_PPARG_SREBF1 = 9.28398;

x_I_PPARG = 10;

% Lambdas
l_HNF4A_HNF4A = 2; %%4.738871;
l_PPARG_PPARG = 9.514515;
l_SREBF1_SREBF1 = 2;%1.302063;  %%#% 2
l_HNF1A_HNF1A = 2; %%4.3344629999999995;

l_HNF4A_HNF1A = 5.328115;
l_HNF1A_HNF4A = 7.737933;
l_SREBF1_PPARG =  6.512859;      %%#% 3
l_PPARG_SREBF1 = 2.7292400000000003;    %%#% 3-6

l_SREBF1_HNF4A = 0.115335;  %%#% 0.4-0.5
l_HNF1A_PPARG = 0.680057;       %%#% 0.5-0.7

l_I_PPARG = 1;

% co-operativity
n_HNF4A_HNF4A = 5.0;
n_HNF4A_HNF1A = 6.0;
n_HNF1A_HNF4A = 3.0;  %%#% 3-4
n_HNF1A_HNF1A = 4.0;

n_HNF1A_PPARG  = 5.0;   %%#% 3-5
n_SREBF1_HNF4A = 2.0;  %%#% 2-4-6

n_PPARG_PPARG  = 5;%2.0;
n_SREBF1_PPARG = 2.0;   %%#% 2-4
n_PPARG_SREBF1 = 2.0;   %%#% 1?? not likely
n_SREBF1_SREBF1= 2.0;   %%#% 2


% all the notation used here is from first to second, that the value
% corresponds to the effect of the first gene on the second



%% Calling of hill function:
H_HNF4A_HNF4A=hill(x(1),x_HNF4A_HNF4A,l_HNF4A_HNF4A,n_HNF4A_HNF4A);
H_HNF1A_HNF4A=hill(x(2),x_HNF1A_HNF4A,l_HNF1A_HNF4A,n_HNF1A_HNF4A);
H_HNF4A_HNF1A=hill(x(1),x_HNF4A_HNF1A,l_HNF4A_HNF1A,n_HNF4A_HNF1A);
H_HNF1A_HNF1A=hill(x(2),x_HNF1A_HNF1A,l_HNF1A_HNF1A,n_HNF1A_HNF1A);

H_SREBF1_HNF4A=hill(x(4),x_SREBF1_HNF4A,l_SREBF1_HNF4A,n_SREBF1_HNF4A);
H_HNF1A_PPARG=hill(x(2),x_HNF1A_PPARG,l_HNF1A_PPARG,n_HNF1A_PPARG);

H_PPARG_PPARG=hill(x(3),x_PPARG_PPARG,l_PPARG_PPARG,n_PPARG_PPARG);
H_SREBF1_PPARG=hill(x(4),x_SREBF1_PPARG,l_SREBF1_PPARG,n_SREBF1_PPARG);
H_PPARG_SREBF1=hill(x(3),x_PPARG_SREBF1,l_PPARG_SREBF1,n_PPARG_SREBF1);
H_SREBF1_SREBF1=hill(x(4),x_SREBF1_SREBF1,l_SREBF1_SREBF1,n_SREBF1_SREBF1);

%% equations:
dxdt(1) = g_HNF4A*H_SREBF1_HNF4A*H_HNF4A_HNF4A*H_HNF1A_HNF4A - k_HNF4A*x(1);
dxdt(2) = g_HNF1A*H_HNF4A_HNF1A*H_HNF1A_HNF1A - k_HNF1A*x(2);
dxdt(3) = g_PPARG*H_HNF1A_PPARG*H_SREBF1_PPARG*H_PPARG_PPARG - k_PPARG*x(3);
dxdt(4) = g_SREBF1*H_PPARG_SREBF1*H_SREBF1_SREBF1 - k_SREBF1*x(4);



% g_A = 55.75;
% g_B = 25.4;
% k_A = bifurcation_parameter;%0.47;
% k_B = 0.21;
% x_A_A = 6.75;
% l_A_A = 10.27;
% n_A_A = 1;
% x_A_B = 8.32;
% l_A_B = 0.011;
% n_A_B = 2;
% x_B_A = 7.84;
% l_B_A = 0.01;
% n_B_A = 3;
% x_B_B = 4.84;
% l_B_B = 28.84;
% n_B_B = 2;
% 
% %% Calling of hill function:
% H_A_A=hill(x(1),x_A_A,l_A_A,n_A_A);
% H_A_B=hill(x(1),x_A_B,l_A_B,n_A_B);
% H_B_A=hill(x(2),x_B_A,l_B_A,n_B_A);
% H_B_B=hill(x(2),x_B_B,l_B_B,n_B_B);
% 
% %% equations:
% dxdt(1) =  g_A*H_A_A*H_B_A - k_A*x(1);
% dxdt(2) = g_B*H_B_B*H_A_B - k_B*x(2);
end
