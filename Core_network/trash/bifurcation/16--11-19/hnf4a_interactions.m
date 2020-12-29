function dxdt = interactions(t,x,bifurcation_parameter)
%% X1: HNF4A, X2: HNF1A, X3:PPARG, X4:SREBF1,
dxdt = zeros(4,1);

%%% production rate
g_HNF4A = 0.4540700422860051; %g/k = 9.2;   2.3-6.9
g_HNF1A = 0.05705569053323611; %g/k = 1.15; 0.29 - 0.87
g_PPARG = 0.09283939658384676; %g/k = 0.56 - 0.238; 
g_SREBF1= 1.1058054056141339; %g/k = 7;     1.75 - 5.25

% Degradation Rate
k_HNF4A = 0.049119799999999997;   %% #% 0.0479
k_HNF1A = 0.049445600000000006; %% #% 0.05
k_PPARG = bifurcation_parameter;%0.07670399999999999;  %% #% 0.06     0.11 to 0.26
k_SREBF1 = 0.15984530000000001; %% #% 0.1732

%Thresholds
x_HNF4A_HNF4A = 3;%7.445485000000001;
x_SREBF1_SREBF1 = 5.25; %4.6490279999999995;      %------------------------------fixed for now
x_PPARG_PPARG = 1;%3;%6.32011;
x_HNF1A_HNF1A = 0.8;     %15.796237;        ------------------------fixed for now

x_SREBF1_HNF4A = 5; %7.067421;                         %-----------------------------fixed for now
x_HNF1A_PPARG = 0.6740759999999998;

x_HNF4A_HNF1A = 5.108258999999999;
x_HNF1A_HNF4A = 0.8; %14.169678;             -----------------------------fixed for now
x_SREBF1_PPARG = 5.25; %2.6013288;               -----------------
x_PPARG_SREBF1 = 3;%9.28398;

x_I_PPARG = 10;

% Lambdas
l_HNF4A_HNF4A = 4; %2; %%4.738871;                       --------------------
l_PPARG_PPARG = 6;%9.514515;
l_SREBF1_SREBF1 = 4; %2;%1.302063;  %%#% 2           --------------------
l_HNF1A_HNF1A = 2; %%4.3344629999999995;

l_HNF4A_HNF1A = 10 ;%5.328115;
l_HNF1A_HNF4A = 4; %7.737933; %%#% 3-4               ---------------------
l_SREBF1_PPARG =  6;% 6.512859;      %%#% 3          ---------------fixed for now
l_PPARG_SREBF1 = 6;%3.7292400000000003;    %%#% 3-6       ------------- fixed for now

l_SREBF1_HNF4A = 0.415335;  %%#% 0.4-0.5             --------------------- fixed for now
l_HNF1A_PPARG = 0.01;%0.680057;       %%#% 0.5-0.7

l_I_PPARG = 1;

% co-operativity
n_HNF4A_HNF4A = 4.0; %5 ------------
n_HNF4A_HNF1A = 6.0; %6 ------------
n_HNF1A_HNF4A = 3.0;
n_HNF1A_HNF1A = 4.0;

n_HNF1A_PPARG  = 4; %5.0;   %%#% 3-5
n_SREBF1_HNF4A = 2.0;  %%#% 2-4-6

n_PPARG_PPARG  = 5.0;   %2.0;
n_SREBF1_PPARG = 2.0;   %%#% 2 -4
n_PPARG_SREBF1 = 2.0;   %%#% 1?? not likely
n_SREBF1_SREBF1= 2.0;   %%#% 2

n_I_PPARG = 2;
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

% equations:
dxdt(2) = g_HNF1A*H_HNF4A_HNF1A*H_HNF1A_HNF1A - k_HNF1A*x(2);
dxdt(4) = g_SREBF1*H_PPARG_SREBF1*H_SREBF1_SREBF1 - k_SREBF1*x(4); 

%% nullcline wrt d(hn4a)dt=0:
dxdt(1) =0;
dxdt(3) = g_PPARG*H_HNF1A_PPARG*H_SREBF1_PPARG*H_PPARG_PPARG - k_PPARG*x(3);
end
