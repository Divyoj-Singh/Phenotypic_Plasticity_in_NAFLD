function [x,v,s,h,f] = interactions
% this function contains all the parameters and calls the function equations.
% it is called by the function bifurcation.    
curdir = pwd;
init;
cd(curdir);

%% Parameters to decide the fineness of the bifurcation plot, step size, maximum number of steps etc.
opt = contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',50000);
opt=contset(opt,'MinStepsize',0.001);
opt=contset(opt,'MaxStepsize',0.1);
opt=contset(opt,'Eigenvalues',1);

%%
ap = 8; %describes the index of parameter for which the bifurcation is drawn using the init_EP_EP function. Currently, ap=8, thus bifurcation parameter is PPARG degradation rate.
handles = feval(@equations);
tspan = 0:1:5000;

% initial condition
x_start = [1 1 1000 1000];
I = 15000;

%%% Production Rates:
g_HNF4A = 0.4540700422860051; 
g_HNF1A = 0.05705569053323611;
g_PPARG = 0.06283939658384676; 
g_SREBF1= 1.1058054056141339; 

% Degradation Rates:
k_HNF4A = 0.049119799999999997;
k_HNF1A = 0.049445600000000006;
k_PPARG = 0.08;
k_SREBF1 = 0.15984530000000001; 

% Thresholds
x_HNF4A_HNF4A = 3;
x_SREBF1_SREBF1 = 5.25;
x_PPARG_PPARG = 6.32011;
x_HNF1A_HNF1A = 0.8;

x_SREBF1_HNF4A = 5; 
x_HNF1A_PPARG = 0.6740759999999998;

x_HNF4A_HNF1A = 5.108258999999999;
x_HNF1A_HNF4A = 0.8;
x_SREBF1_PPARG = 5.25;
x_PPARG_SREBF1 = 9.28398;

x_I_PPARG = 10;

% Lambdas (fold change)
l_HNF4A_HNF4A = 4;
l_PPARG_PPARG = 9.514515;
l_SREBF1_SREBF1 = 4;
l_HNF1A_HNF1A = 2;

l_HNF4A_HNF1A = 5.328115;
l_HNF1A_HNF4A = 4; 
l_SREBF1_PPARG =  3;
l_PPARG_SREBF1 = 3.7292400000000003; 

l_SREBF1_HNF4A = 0.415335;
l_HNF1A_PPARG = 0.680057;

l_I_PPARG = 1;

% co-operativity
n_HNF4A_HNF4A = 4.0;
n_HNF4A_HNF1A = 4.0;
n_HNF1A_HNF4A = 3.0;
n_HNF1A_HNF1A = 4.0;

n_HNF1A_PPARG  = 4;
n_SREBF1_HNF4A = 2.0;

n_PPARG_PPARG  = 5.0; 
n_SREBF1_PPARG = 2.0;
n_PPARG_SREBF1 = 2.0;
n_SREBF1_SREBF1= 2.0;
n_I_PPARG = 2;

% Calculating steady state for given initial condition: 
[t,x_time] = ode15s(@(t,protein_array)handles{2}(t,protein_array,I,g_HNF4A,g_HNF1A,g_PPARG,g_SREBF1,k_HNF4A,k_HNF1A ,k_PPARG ,k_SREBF1 ,n_HNF4A_HNF4A ,n_HNF4A_HNF1A ,n_HNF1A_HNF4A ,n_HNF1A_HNF1A,n_HNF1A_PPARG,n_SREBF1_HNF4A,n_PPARG_PPARG,n_SREBF1_PPARG,n_PPARG_SREBF1,n_SREBF1_SREBF1,l_HNF4A_HNF4A ,l_HNF4A_HNF1A ,l_HNF1A_HNF4A ,l_HNF1A_HNF1A ,l_SREBF1_HNF4A,l_HNF1A_PPARG,l_PPARG_PPARG,l_SREBF1_PPARG,l_PPARG_SREBF1,l_SREBF1_SREBF1,x_HNF4A_HNF4A,x_HNF4A_HNF1A ,x_HNF1A_HNF4A ,x_HNF1A_HNF1A ,x_SREBF1_HNF4A,x_HNF1A_PPARG,x_PPARG_PPARG,x_SREBF1_PPARG,x_PPARG_SREBF1,x_SREBF1_SREBF1,x_I_PPARG,l_I_PPARG,n_I_PPARG),tspan,x_start);
x_init = x_time(end,:)'

% Drawing bifurcation using a continuation method:
[x0,v0] = init_EP_EP(@equations,x_init,[I g_HNF4A g_HNF1A g_PPARG g_SREBF1 k_HNF4A k_HNF1A k_PPARG k_SREBF1 n_HNF4A_HNF4A n_HNF4A_HNF1A n_HNF1A_HNF4A n_HNF1A_HNF1A n_HNF1A_PPARG n_SREBF1_HNF4A n_PPARG_PPARG n_SREBF1_PPARG n_PPARG_SREBF1 n_SREBF1_SREBF1 l_HNF4A_HNF4A l_HNF4A_HNF1A l_HNF1A_HNF4A l_HNF1A_HNF1A l_SREBF1_HNF4A l_HNF1A_PPARG l_PPARG_PPARG l_SREBF1_PPARG l_PPARG_SREBF1 l_SREBF1_SREBF1 x_HNF4A_HNF4A x_HNF4A_HNF1A x_HNF1A_HNF4A x_HNF1A_HNF1A x_SREBF1_HNF4A x_HNF1A_PPARG x_PPARG_PPARG x_SREBF1_PPARG x_PPARG_SREBF1 x_SREBF1_SREBF1,x_I_PPARG,l_I_PPARG,n_I_PPARG],ap);%%%% parameters(n) needs to be an 1d array
[x,v,s,h,f] = cont(@equilibrium, x0, v0,opt);

end
