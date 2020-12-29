function out = equations_hnf4a
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,protein_array,I,g_HNF4A,g_HNF1A,g_PPARG,g_SREBF1,k_HNF4A,k_HNF1A ,k_PPARG ,k_SREBF1 ,n_HNF4A_HNF4A ,n_HNF4A_HNF1A ,n_HNF1A_HNF4A ,n_HNF1A_HNF1A,n_HNF1A_PPARG,n_SREBF1_HNF4A,n_PPARG_PPARG,n_SREBF1_PPARG,n_PPARG_SREBF1,n_SREBF1_SREBF1,l_HNF4A_HNF4A ,l_HNF4A_HNF1A ,l_HNF1A_HNF4A ,l_HNF1A_HNF1A ,l_SREBF1_HNF4A,l_HNF1A_PPARG,l_PPARG_PPARG,l_SREBF1_PPARG,l_PPARG_SREBF1,l_SREBF1_SREBF1,x_HNF4A_HNF4A,x_HNF4A_HNF1A ,x_HNF1A_HNF4A ,x_HNF1A_HNF1A ,x_SREBF1_HNF4A,x_HNF1A_PPARG,x_PPARG_PPARG,x_SREBF1_PPARG,x_PPARG_SREBF1,x_SREBF1_SREBF1,x_I_PPARG,l_I_PPARG,n_I_PPARG)


%% Calling of hill function:
H_HNF4A_HNF4A=hill(I,x_HNF4A_HNF4A,l_HNF4A_HNF4A,n_HNF4A_HNF4A);
H_HNF1A_HNF4A=hill(protein_array(1),x_HNF1A_HNF4A,l_HNF1A_HNF4A,n_HNF1A_HNF4A);
H_HNF4A_HNF1A=hill(I,x_HNF4A_HNF1A,l_HNF4A_HNF1A,n_HNF4A_HNF1A);
H_HNF1A_HNF1A=hill(protein_array(1),x_HNF1A_HNF1A,l_HNF1A_HNF1A,n_HNF1A_HNF1A);

H_SREBF1_HNF4A=hill(protein_array(3),x_SREBF1_HNF4A,l_SREBF1_HNF4A,n_SREBF1_HNF4A);
H_HNF1A_PPARG=hill(protein_array(1),x_HNF1A_PPARG,l_HNF1A_PPARG,n_HNF1A_PPARG);

H_PPARG_PPARG=hill(protein_array(2),x_PPARG_PPARG,l_PPARG_PPARG,n_PPARG_PPARG);
H_SREBF1_PPARG=hill(protein_array(3),x_SREBF1_PPARG,l_SREBF1_PPARG,n_SREBF1_PPARG);
H_PPARG_SREBF1=hill(protein_array(2),x_PPARG_SREBF1,l_PPARG_SREBF1,n_PPARG_SREBF1);
H_SREBF1_SREBF1=hill(protein_array(3),x_SREBF1_SREBF1,l_SREBF1_SREBF1,n_SREBF1_SREBF1);


%% equations:
% dxdt(1) = g_HNF4A*H_SREBF1_HNF4A*H_HNF4A_HNF4A*H_HNF1A_HNF4A - k_HNF4A*x(1);
% dxdt(2) = g_HNF1A*H_HNF4A_HNF1A*H_HNF1A_HNF1A - k_HNF1A*x(2);
% dxdt(3) = g_PPARG*H_HNF1A_PPARG*H_SREBF1_PPARG*H_PPARG_PPARG - k_PPARG*x(3);
% dxdt(4) = g_SREBF1*H_PPARG_SREBF1*H_SREBF1_SREBF1 - k_SREBF1*x(4);

dydt=[ g_HNF1A*H_HNF4A_HNF1A*H_HNF1A_HNF1A - k_HNF1A*protein_array(1);
    g_PPARG*H_HNF1A_PPARG*H_SREBF1_PPARG*H_PPARG_PPARG - k_PPARG*protein_array(2);
    g_SREBF1*H_PPARG_SREBF1*H_SREBF1_SREBF1 - k_SREBF1*protein_array(3);
];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(equations_hnf4a);
y0=[0,0,0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,kmgr,kgr,gmgr,ggr,gr0mz,z0mgr,ngrmz,nzmgr,lamdagrmz,lamdazmgr)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,kmgr,kgr,gmgr,ggr,gr0mz,z0mgr,ngrmz,nzmgr,lamdagrmz,lamdazmgr)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,kmgr,kgr,gmgr,ggr,gr0mz,z0mgr,ngrmz,nzmgr,lamdagrmz,lamdazmgr)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,kmgr,kgr,gmgr,ggr,gr0mz,z0mgr,ngrmz,nzmgr,lamdagrmz,lamdazmgr)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,kmgr,kgr,gmgr,ggr,gr0mz,z0mgr,ngrmz,nzmgr,lamdagrmz,lamdazmgr)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,kmgr,kgr,gmgr,ggr,gr0mz,z0mgr,ngrmz,nzmgr,lamdagrmz,lamdazmgr)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,kmgr,kgr,gmgr,ggr,gr0mz,z0mgr,ngrmz,nzmgr,lamdagrmz,lamdazmgr)
