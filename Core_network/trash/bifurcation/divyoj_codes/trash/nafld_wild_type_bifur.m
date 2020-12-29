function [x,v,s,h,f] = nafld_wild_type_bifur

curdir = pwd;
init;
cd(curdir);

opt = contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',50000);
opt=contset(opt,'MinStepsize',0.001);
opt=contset(opt,'MaxStepsize',10);
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
ap = 1; %describes the index of parameter for which the bifurcation is drawn using the init_EP_EP function. Currently, ap=1, thus bifurcation parameter is s (SNAIL levels)
handles = feval(@nafld_wild_type);
tspan = 0:1:5000;

% initial condition
x_start = [1 1 10000 10  000];
I = 0;

%%% production rate
g_HNF4A= 12.073225;%42.420912
g_HNF1A= 8.29212375; %30.865732;
g_PPARG= 6.7606944; %56.121306;
g_SREBF1= 23.57029275; %94.768549;

% Degradation Rate
k_HNF4A = 0.0479; %0.178136
k_HNF1A = 0.05; %0.208767
k_PPARG = 0.31164; %0.135633   for bistable, k_PPARG=0.05;
k_SREBF1 = 0.1732; %0.377047

%Thresholds
x_HNF4A_HNF4A = 5; %4.992633;  %g/k=400;
x_SREBF1_SREBF1= 15; %13.425561;  %g/k=335;
x_PPARG_PPARG= 1; %0.774717;  %g/k=600;
x_HNF1A_HNF1A = 5; %6.060171;   %g/k=400;

x_SREBF1_HNF4A= 200; %18.985005; g/k=335;
x_HNF1A_PPARG= 80; %22.45355;   g/k=400;

x_HNF4A_HNF1A = 400; %6.152792;   g/k=400;
x_HNF1A_HNF4A = 80; %18.407056;   g/k=400;
x_SREBF1_PPARG= 500; %9.196498;    g/k=335;
x_PPARG_SREBF1= 300; %4.629643;     g/k=600;

x_I_PPARG = 20;

% Lambdas
l_HNF4A_HNF4A = 1; %99.909396;
l_PPARG_PPARG = 1; %55.233684;
l_SREBF1_SREBF1 = 1; %4.043209;
l_HNF1A_HNF1A = 1; %23.675277;

l_HNF4A_HNF1A = 4; %50.888983;
l_HNF1A_HNF4A = 4; %9.393584;
l_SREBF1_PPARG = 4; %20.506591;
l_PPARG_SREBF1 = 4; %88.014212;

l_SREBF1_HNF4A = 0.1; %0.019458;
l_HNF1A_PPARG = 0.1; %0.016129;

l_I_PPARG = 2;

% co-operativity
n_HNF4A_HNF4A = 3;
n_HNF4A_HNF1A = 4;
n_HNF1A_HNF4A = 5;
n_HNF1A_HNF1A = 3;

n_HNF1A_PPARG  = 4;
n_SREBF1_HNF4A = 4;

n_PPARG_PPARG  = 2;
n_SREBF1_PPARG = 5;
n_PPARG_SREBF1 = 4;
n_SREBF1_SREBF1= 5; % 5;

n_I_PPARG = 2;

% %% Paramaters
% % n is the id of parameter set.
% 
% % production rate
% g_HNF4A=parameters_all.Prod_of_HNF4A(n);
% g_HNF1A=parameters_all.Prod_of_HNF1A(n);
% g_PPARG=parameters_all.Prod_of_PPARG(n);
% 
% g_SREBF1=parameters_all.Prod_of_SREBF1(n);
% 
% % Degradation Rate
% k_HNF4A =parameters_all.Deg_of_HNF4A(n);
% k_HNF1A = parameters_all.Deg_of_HNF1A(n);
% k_PPARG = parameters_all.Deg_of_PPARG(n);
% k_PPARG = 0.335633;
% k_SREBF1 = parameters_all.Deg_of_SREBF1(n);
% 
% % co-operativity
% n_HNF4A_HNF4A = parameters_all.Num_of_HNF4AToHNF4A(n);
% n_HNF4A_HNF1A =parameters_all.Num_of_HNF4AToHNF1A(n);
% n_HNF1A_HNF4A =parameters_all.Num_of_HNF1AToHNF4A(n);
% n_HNF1A_HNF1A=parameters_all.Num_of_HNF1AToHNF1A(n);
% 
% n_HNF1A_PPARG=parameters_all.Num_of_HNF1AToPPARG(n);
% n_SREBF1_HNF4A=parameters_all.Num_of_SREBF1ToHNF4A(n);
% 
% n_PPARG_PPARG=parameters_all.Num_of_PPARGToPPARG(n);
% n_SREBF1_PPARG=parameters_all.Num_of_SREBF1ToPPARG(n);
% n_PPARG_SREBF1=parameters_all.Num_of_PPARGToSREBF1(n);
% n_SREBF1_SREBF1=parameters_all.Num_of_SREBF1ToSREBF1(n);
% 
% 
% % Lambdas
% l_HNF4A_HNF4A = parameters_all.Act_of_HNF4AToHNF4A(n);
% l_HNF4A_HNF1A =parameters_all.Act_of_HNF4AToHNF1A(n);
% l_HNF1A_HNF4A =parameters_all.Act_of_HNF1AToHNF4A(n) ;
% l_HNF1A_HNF1A = parameters_all.Act_of_HNF1AToHNF1A(n);
% 
% l_SREBF1_HNF4A=parameters_all.Inh_of_SREBF1ToHNF4A(n);
% l_HNF1A_PPARG=parameters_all.Inh_of_HNF1AToPPARG(n);
% 
% l_PPARG_PPARG=parameters_all.Act_of_PPARGToPPARG(n);
% l_SREBF1_PPARG=parameters_all.Act_of_SREBF1ToPPARG(n);
% l_PPARG_SREBF1=parameters_all.Act_of_PPARGToSREBF1(n);
% l_SREBF1_SREBF1=parameters_all.Act_of_SREBF1ToSREBF1(n);
% 
% x_HNF4A_HNF4A = parameters_all.Trd_of_HNF4AToHNF4A(n);
% x_HNF4A_HNF1A =parameters_all.Trd_of_HNF4AToHNF1A(n);
% x_HNF1A_HNF4A =parameters_all.Trd_of_HNF1AToHNF4A(n); 
% x_HNF1A_HNF1A = parameters_all.Trd_of_HNF1AToHNF1A(n);
% 
% x_SREBF1_HNF4A=parameters_all.Trd_of_SREBF1ToHNF4A(n);
% x_HNF1A_PPARG=parameters_all.Trd_of_HNF1AToPPARG(n);
% 
% x_PPARG_PPARG=parameters_all.Trd_of_PPARGToPPARG(n);
% x_SREBF1_PPARG=parameters_all.Trd_of_SREBF1ToPPARG(n);
% x_PPARG_SREBF1=parameters_all.Trd_of_PPARGToSREBF1(n);
% x_SREBF1_SREBF1=parameters_all.Trd_of_SREBF1ToSREBF1(n);
%%

%calculating steady state for given initial condition 
[t,x_time] = ode15s(@(t,protein_array)handles{2}(t,protein_array,I,g_HNF4A,g_HNF1A,g_PPARG,g_SREBF1,k_HNF4A,k_HNF1A ,k_PPARG ,k_SREBF1 ,n_HNF4A_HNF4A ,n_HNF4A_HNF1A ,n_HNF1A_HNF4A ,n_HNF1A_HNF1A,n_HNF1A_PPARG,n_SREBF1_HNF4A,n_PPARG_PPARG,n_SREBF1_PPARG,n_PPARG_SREBF1,n_SREBF1_SREBF1,l_HNF4A_HNF4A ,l_HNF4A_HNF1A ,l_HNF1A_HNF4A ,l_HNF1A_HNF1A ,l_SREBF1_HNF4A,l_HNF1A_PPARG,l_PPARG_PPARG,l_SREBF1_PPARG,l_PPARG_SREBF1,l_SREBF1_SREBF1,x_HNF4A_HNF4A,x_HNF4A_HNF1A ,x_HNF1A_HNF4A ,x_HNF1A_HNF1A ,x_SREBF1_HNF4A,x_HNF1A_PPARG,x_PPARG_PPARG,x_SREBF1_PPARG,x_PPARG_SREBF1,x_SREBF1_SREBF1,x_I_PPARG,l_I_PPARG,n_I_PPARG),tspan,x_start);
x_init = x_time(end,:)'

%drawing bifurcation using a continuation method
[x0,v0] = init_EP_EP(@nafld_wild_type,x_init,[I g_HNF4A g_HNF1A g_PPARG g_SREBF1 k_HNF4A k_HNF1A k_PPARG k_SREBF1 n_HNF4A_HNF4A n_HNF4A_HNF1A n_HNF1A_HNF4A n_HNF1A_HNF1A n_HNF1A_PPARG n_SREBF1_HNF4A n_PPARG_PPARG n_SREBF1_PPARG n_PPARG_SREBF1 n_SREBF1_SREBF1 l_HNF4A_HNF4A l_HNF4A_HNF1A l_HNF1A_HNF4A l_HNF1A_HNF1A l_SREBF1_HNF4A l_HNF1A_PPARG l_PPARG_PPARG l_SREBF1_PPARG l_PPARG_SREBF1 l_SREBF1_SREBF1 x_HNF4A_HNF4A x_HNF4A_HNF1A x_HNF1A_HNF4A x_HNF1A_HNF1A x_SREBF1_HNF4A x_HNF1A_PPARG x_PPARG_PPARG x_SREBF1_PPARG x_PPARG_SREBF1 x_SREBF1_SREBF1,x_I_PPARG,l_I_PPARG,n_I_PPARG],ap);%%%% parameters(n) needs to be an 1d array
[x,v,s,h,f] = cont(@equilibrium, x0, v0,opt);

end


