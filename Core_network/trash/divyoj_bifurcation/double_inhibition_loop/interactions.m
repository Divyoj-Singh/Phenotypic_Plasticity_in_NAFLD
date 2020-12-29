function dxdt = interactions(t,x,bifurcation_parameter)
%% X1: HNF4A, X2: HNF1A, X3:PPARG, X4:SREBF1,
dxdt = zeros(2,1);

g_A = 55.75;
g_B = 25.4;
k_A = bifurcation_parameter;%0.47;
k_B = 0.21;
x_A_A = 6.75;
l_A_A = 10.27;
n_A_A = 1;
x_A_B = 8.32;
l_A_B = 0.011;
n_A_B = 2;
x_B_A = 7.84;
l_B_A = 0.01;
n_B_A = 3;
x_B_B = 4.84;
l_B_B = 28.84;
n_B_B = 2;

%% Calling of hill function:
H_A_A=hill(x(1),x_A_A,l_A_A,n_A_A);
H_A_B=hill(x(1),x_A_B,l_A_B,n_A_B);
H_B_A=hill(x(2),x_B_A,l_B_A,n_B_A);
H_B_B=hill(x(2),x_B_B,l_B_B,n_B_B);

%% equations:
dxdt(1) =  g_A*H_A_A*H_B_A - k_A*x(1);
dxdt(2) = g_B*H_B_B*H_A_B - k_B*x(2);
end
