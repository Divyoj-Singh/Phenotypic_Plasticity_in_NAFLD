function dxdt = interactions(t,x,s,n)
%% X1: HNF4A, X2: HNF1A, X3:PPARG, X4:SREBF1,
dxdt = zeros(4,1);

% all the notation used here is from first to second, that the value
% corresponds to the effect of the first gene on the second

%% Paramaters
% n is the id of parameter set.

% production rate
g_HNF4A=s.Prod_of_HNF4A(n);
g_HNF1A=s.Prod_of_HNF1A(n);
g_PPARG=s.Prod_of_PPARG(n);
g_SREBF1=s.Prod_of_SREBF1(n);

% Degradation Rate
k_HNF4A =s.Deg_of_HNF4A(n);
k_HNF1A = s.Deg_of_HNF1A(n);
k_PPARG = s.Deg_of_PPARG(n);
k_SREBF1 = s.Deg_of_SREBF1(n);

% co-operativity
n_HNF4A_HNF4A = s.Num_of_HNF4AToHNF4A(n);
n_HNF4A_HNF1A =s.Num_of_HNF4AToHNF1A(n);
n_HNF1A_HNF4A =s.Num_of_HNF1AToHNF4A(n);
n_HNF1A_HNF1A=s.Num_of_HNF1AToHNF1A(n);

n_HNF1A_PPARG=s.Num_of_HNF1AToPPARG(n);
n_SREBF1_HNF4A=s.Num_of_SREBF1ToHNF4A(n);

n_PPARG_PPARG=s.Num_of_PPARGToPPARG(n);
n_SREBF1_PPARG=s.Num_of_SREBF1ToPPARG(n);
n_PPARG_SREBF1=s.Num_of_PPARGToSREBF1(n);
n_SREBF1_SREBF1=s.Num_of_SREBF1ToSREBF1(n);


% Lambdas(fold change)
l_HNF4A_HNF4A = s.Act_of_HNF4AToHNF4A(n);
l_HNF4A_HNF1A =s.Act_of_HNF4AToHNF1A(n);
l_HNF1A_HNF4A =s.Act_of_HNF1AToHNF4A(n) ;
l_HNF1A_HNF1A = s.Act_of_HNF1AToHNF1A(n);

l_SREBF1_HNF4A=s.Inh_of_SREBF1ToHNF4A(n);
l_HNF1A_PPARG=s.Inh_of_HNF1AToPPARG(n);

l_PPARG_PPARG=s.Act_of_PPARGToPPARG(n);
l_SREBF1_PPARG=s.Act_of_SREBF1ToPPARG(n);
l_PPARG_SREBF1=s.Act_of_PPARGToSREBF1(n);
l_SREBF1_SREBF1=s.Act_of_SREBF1ToSREBF1(n);

% constant in hill function:

x_HNF4A_HNF4A = s.Trd_of_HNF4AToHNF4A(n);
x_HNF4A_HNF1A =s.Trd_of_HNF4AToHNF1A(n);
x_HNF1A_HNF4A =s.Trd_of_HNF1AToHNF4A(n); 
x_HNF1A_HNF1A = s.Trd_of_HNF1AToHNF1A(n);

x_SREBF1_HNF4A=s.Trd_of_SREBF1ToHNF4A(n);
x_HNF1A_PPARG=s.Trd_of_HNF1AToPPARG(n);

x_PPARG_PPARG=s.Trd_of_PPARGToPPARG(n);
x_SREBF1_PPARG=s.Trd_of_SREBF1ToPPARG(n);
x_PPARG_SREBF1=s.Trd_of_PPARGToSREBF1(n);
x_SREBF1_SREBF1=s.Trd_of_SREBF1ToSREBF1(n);


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
end
