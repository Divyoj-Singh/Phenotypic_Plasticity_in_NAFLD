function [x,v,s,h,f] = nafld_wild_type_bifur

curdir = pwd;
init;
cd(curdir);

opt = contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',50000);
opt=contset(opt,'MinStepsize',0.001);
opt=contset(opt,'MaxStepsize',0.1);
opt=contset(opt,'Eigenvalues',1);

% %% reading the parameters file:
% % change the parameter file path here:
% parameters_all= tdfread('./../state_2_HL_LH_parameter_subset.txt','\t');
% n = 1;
% parameter=[];
% fn = fieldnames(parameters_all);
% for k=3:numel(fn)
%     if( isnumeric(parameters_all.(fn{k})) )
%         parameter=[parameter,parameters_all.(fn{k})(n)];
%     end
% end
 
%%
ap = 8; %describes the index of parameter for which the bifurcation is drawn using the init_EP_EP function. Currently, ap=1, thus bifurcation parameter is s (SNAIL levels)
handles = feval(@nafld_wild_type);
tspan = 0:1:5000;

% initial condition
x_start = [1 1 1000 1000];
I = 15000;

%%% production rate
g_HNF4A = 0.4540700422860051; %g/k = 9.2;   2.3-6.9
g_HNF1A = 0.05705569053323611;%g/k = 1.15; 0.29 - 0.87
g_PPARG = 0.09583939658384676; %g/k = 0.56 - 0.238; %% x0 0.23-0.71
g_SREBF1= 1.1058054056141339; %g/k = 7;     1.75 - 5.25

% Degradation Rate
k_HNF4A = 0.049119799999999997;   %% #% 0.0479
k_HNF1A = 0.049445600000000006; %% #% 0.05
k_PPARG = 0.001;%0.07670399999999999;  %% #% 0.06     0.11 to 0.26  according to 3hrs to 12 hrs -- 0.23-0.0575 
k_SREBF1 = 0.15984530000000001; %% #% 0.1732

%Thresholds
x_HNF4A_HNF4A = 2.3;%7.445485000000001;
x_SREBF1_SREBF1 = 5.25; %4.6490279999999995;      %------------------------------fixed for now
x_PPARG_PPARG = 0.8;%6.32011;
x_HNF1A_HNF1A = 0.5;     %15.796237;        ------------------------fixed for now

x_SREBF1_HNF4A = 5; %7.067421;                         %-----------------------------fixed for now
x_HNF1A_PPARG = 0.6740759999999998;

x_HNF4A_HNF1A = 7;%5.108258999999999;
x_HNF1A_HNF4A = 0.8; %14.169678;             -----------------------------fixed for now
x_SREBF1_PPARG = 5.25; %2.6013288;               -----------------
x_PPARG_SREBF1 = 1;%%%9.28398;

x_I_PPARG = 10;

% Lambdas
l_HNF4A_HNF4A = 5; %2; %%4.738871;                       --------------------
l_PPARG_PPARG = 5;%%9.514515;
l_SREBF1_SREBF1 = 2; %2;%1.302063;  %%#% 2           --------------------
l_HNF1A_HNF1A = 7; %%4.3344629999999995;

l_HNF4A_HNF1A = 9;%5.328115;
l_HNF1A_HNF4A = 1.5;%7.737933; %%#% 3-4               ---------------------
l_SREBF1_PPARG =  3;% 6.512859;      %%#% 3          ---------------fixed for now
l_PPARG_SREBF1 = 3;%3.7292400000000003;    %%#% 3-6       ------------- fixed for now

l_SREBF1_HNF4A = 0.10335;  %%#% 0.4-0.5             --------------------- fixed for now
l_HNF1A_PPARG = 0.6800572;       %%#% 0.5-0.7

l_I_PPARG = 1;

% co-operativity
n_HNF4A_HNF4A = 5;%4.0; %5 ------------
n_HNF4A_HNF1A = 4.0; %6 ------------
n_HNF1A_HNF4A = 4;%3.0;
n_HNF1A_HNF1A = 4.0;

n_HNF1A_PPARG  = 5; %5.0;   %%#% 3-5
n_SREBF1_HNF4A = 2.0;  %%#% 2-4-6

n_PPARG_PPARG  = 4.0;   %2.0;
n_SREBF1_PPARG = 2.0;   %--2.0;   %%#% 2 -4
n_PPARG_SREBF1 = 2.0;   %%#% 1?? not likely
n_SREBF1_SREBF1= 2.0;   %%#% 2

n_I_PPARG = 2;


%calculating steady state for given initial condition 
[t,x_time] = ode15s(@(t,protein_array)handles{2}(t,protein_array,I,g_HNF4A,g_HNF1A,g_PPARG,g_SREBF1,k_HNF4A,k_HNF1A ,k_PPARG ,k_SREBF1 ,n_HNF4A_HNF4A ,n_HNF4A_HNF1A ,n_HNF1A_HNF4A ,n_HNF1A_HNF1A,n_HNF1A_PPARG,n_SREBF1_HNF4A,n_PPARG_PPARG,n_SREBF1_PPARG,n_PPARG_SREBF1,n_SREBF1_SREBF1,l_HNF4A_HNF4A ,l_HNF4A_HNF1A ,l_HNF1A_HNF4A ,l_HNF1A_HNF1A ,l_SREBF1_HNF4A,l_HNF1A_PPARG,l_PPARG_PPARG,l_SREBF1_PPARG,l_PPARG_SREBF1,l_SREBF1_SREBF1,x_HNF4A_HNF4A,x_HNF4A_HNF1A ,x_HNF1A_HNF4A ,x_HNF1A_HNF1A ,x_SREBF1_HNF4A,x_HNF1A_PPARG,x_PPARG_PPARG,x_SREBF1_PPARG,x_PPARG_SREBF1,x_SREBF1_SREBF1,x_I_PPARG,l_I_PPARG,n_I_PPARG),tspan,x_start);
x_init = x_time(end,:)'

%drawing bifurcation using a continuation method
[x0,v0] = init_EP_EP(@nafld_wild_type,x_init,[I g_HNF4A g_HNF1A g_PPARG g_SREBF1 k_HNF4A k_HNF1A k_PPARG k_SREBF1 n_HNF4A_HNF4A n_HNF4A_HNF1A n_HNF1A_HNF4A n_HNF1A_HNF1A n_HNF1A_PPARG n_SREBF1_HNF4A n_PPARG_PPARG n_SREBF1_PPARG n_PPARG_SREBF1 n_SREBF1_SREBF1 l_HNF4A_HNF4A l_HNF4A_HNF1A l_HNF1A_HNF4A l_HNF1A_HNF1A l_SREBF1_HNF4A l_HNF1A_PPARG l_PPARG_PPARG l_SREBF1_PPARG l_PPARG_SREBF1 l_SREBF1_SREBF1 x_HNF4A_HNF4A x_HNF4A_HNF1A x_HNF1A_HNF4A x_HNF1A_HNF1A x_SREBF1_HNF4A x_HNF1A_PPARG x_PPARG_PPARG x_SREBF1_PPARG x_PPARG_SREBF1 x_SREBF1_SREBF1,x_I_PPARG,l_I_PPARG,n_I_PPARG],ap);%%%% parameters(n) needs to be an 1d array
[x,v,s,h,f] = cont(@equilibrium, x0, v0,opt);

end


