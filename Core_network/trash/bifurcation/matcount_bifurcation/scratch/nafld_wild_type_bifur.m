function [x,v,s,h,f] = nafld_wild_type_bifur

curdir = pwd;
init;
cd(curdir);

opt = contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',20000);
opt=contset(opt,'MinStepsize',0.0001);
opt=contset(opt,'MaxStepsize',10);
opt=contset(opt,'Eigenvalues',1);

% % reading the parameters file:
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
g_HNF4A = 20;
g_HNF1A = 20;
g_PPARG = 20;
g_SREBF1= 20;

% Degradation Rate
k_HNF4A = 0.05;
k_HNF1A = 0.05;
k_PPARG = 1;%0.320929;
k_SREBF1 = 0.5;

%Thresholds
x_HNF4A_HNF4A = 200;
x_SREBF1_SREBF1 = 200;
x_PPARG_PPARG = 200;
x_HNF1A_HNF1A = 200;

x_SREBF1_HNF4A = 200;
x_HNF1A_PPARG = 20;

x_HNF4A_HNF1A = 200;
x_HNF1A_HNF4A = 200;
x_SREBF1_PPARG = 200;
x_PPARG_SREBF1 = 200;

x_I_PPARG = 10;

% Lambdas
l_HNF4A_HNF4A = 1;
l_PPARG_PPARG = 1;
l_SREBF1_SREBF1 = 1;
l_HNF1A_HNF1A = 1;

l_SREBF1_HNF4A = 5;
l_HNF1A_PPARG = 5;

l_HNF4A_HNF1A = 1;
l_HNF1A_HNF4A = 1;
l_SREBF1_PPARG = 1;
l_PPARG_SREBF1 = 1;

l_I_PPARG = 1;

% co-operativity
n_HNF4A_HNF4A = 3.0;
n_HNF4A_HNF1A = 3.0;
n_HNF1A_HNF4A = 3.0;
n_HNF1A_HNF1A = 3.0;

n_HNF1A_PPARG  = 3.0;
n_SREBF1_HNF4A = 3.0;

n_PPARG_PPARG  = 3.0;
n_SREBF1_PPARG = 3.0;
n_PPARG_SREBF1 = 3.0;
n_SREBF1_SREBF1= 3.0;

n_I_PPARG = 2;

%calculating steady state for given initial condition 
[t,x_time] = ode15s(@(t,protein_array)handles{2}(t,protein_array,I,g_HNF4A,g_HNF1A,g_PPARG,g_SREBF1,k_HNF4A,k_HNF1A ,k_PPARG ,k_SREBF1 ,n_HNF4A_HNF4A ,n_HNF4A_HNF1A ,n_HNF1A_HNF4A ,n_HNF1A_HNF1A,n_HNF1A_PPARG,n_SREBF1_HNF4A,n_PPARG_PPARG,n_SREBF1_PPARG,n_PPARG_SREBF1,n_SREBF1_SREBF1,l_HNF4A_HNF4A ,l_HNF4A_HNF1A ,l_HNF1A_HNF4A ,l_HNF1A_HNF1A ,l_SREBF1_HNF4A,l_HNF1A_PPARG,l_PPARG_PPARG,l_SREBF1_PPARG,l_PPARG_SREBF1,l_SREBF1_SREBF1,x_HNF4A_HNF4A,x_HNF4A_HNF1A ,x_HNF1A_HNF4A ,x_HNF1A_HNF1A ,x_SREBF1_HNF4A,x_HNF1A_PPARG,x_PPARG_PPARG,x_SREBF1_PPARG,x_PPARG_SREBF1,x_SREBF1_SREBF1,x_I_PPARG,l_I_PPARG,n_I_PPARG),tspan,x_start);
x_init = x_time(end,:)'

%drawing bifurcation using a continuation method
[x0,v0] = init_EP_EP(@nafld_wild_type,x_init,[I g_HNF4A g_HNF1A g_PPARG g_SREBF1 k_HNF4A k_HNF1A k_PPARG k_SREBF1 n_HNF4A_HNF4A n_HNF4A_HNF1A n_HNF1A_HNF4A n_HNF1A_HNF1A n_HNF1A_PPARG n_SREBF1_HNF4A n_PPARG_PPARG n_SREBF1_PPARG n_PPARG_SREBF1 n_SREBF1_SREBF1 l_HNF4A_HNF4A l_HNF4A_HNF1A l_HNF1A_HNF4A l_HNF1A_HNF1A l_SREBF1_HNF4A l_HNF1A_PPARG l_PPARG_PPARG l_SREBF1_PPARG l_PPARG_SREBF1 l_SREBF1_SREBF1 x_HNF4A_HNF4A x_HNF4A_HNF1A x_HNF1A_HNF4A x_HNF1A_HNF1A x_SREBF1_HNF4A x_HNF1A_PPARG x_PPARG_PPARG x_SREBF1_PPARG x_PPARG_SREBF1 x_SREBF1_SREBF1,x_I_PPARG,l_I_PPARG,n_I_PPARG],ap);%%%% parameters(n) needs to be an 1d array
[x,v,s,h,f] = cont(@equilibrium, x0, v0,opt);
end


