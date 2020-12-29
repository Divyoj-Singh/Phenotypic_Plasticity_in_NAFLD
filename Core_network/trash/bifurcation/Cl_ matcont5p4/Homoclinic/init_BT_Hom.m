function [x,v] = init_BT_Hom(odefile, x, s, p, ap, ntst, ncol, eps, extravec)
global homds cds
%
% 1. CHECCK INPUT
%--------------------------------------------
n_par = size(ap,2);
if n_par ~= 2
    error('2 free system parameters are needed');
end

if extravec==zeros(1,3)
    error('at least 1 free homoclinic parameter is needed');
end

if extravec==ones(1,3)
    error('at most 2 free homoclinic parameters are needed');
end

for i=1:length(s)
    if ~isempty(s.data) && isfield(s(i).data,'evec')
        nph = size(s(i).data.evec,2);
    elseif isfield(s(i).data,'x')
        nph = size(s(i).data.x,2);
    else
        
    end
end

if isempty(cds) || ~isfield(cds,'options')
    cds.options = contset();
end
cds.curve = @homoclinic;
curvehandles = feval(cds.curve);
cds.curve_func = curvehandles{1};
cds.curve_jacobian = curvehandles{4};
cds.curve_hessians = curvehandles{5};
homds = [];
%
% 2. INITIALIZE homds
%--------------------------------------------
kfactor = 1.02; % 
eps0=1e-2;
eps1=1e-2;
T = 100;

init_homds(odefile,x,p,ap,ntst,ncol,extravec,T,eps0,eps1,nph);

func_handles = feval(odefile);
symord = 0; 
symordp = 0;

if     ~isempty(func_handles{9}),   symord = 5; 
elseif ~isempty(func_handles{8}),   symord = 4; 
elseif ~isempty(func_handles{7}),   symord = 3; 
elseif ~isempty(func_handles{5}),   symord = 2; 
elseif ~isempty(func_handles{3}),   symord = 1; 
end
if     ~isempty(func_handles{6}),   symordp = 2; 
elseif ~isempty(func_handles{4}),   symordp = 1; 
end
cds.options = contset(cds.options, 'SymDerivative', symord);
cds.options = contset(cds.options, 'SymDerivativeP', symordp);
cds.symjac = 1;
cds.symhess = 0;

homds.odefile = odefile;
homds.func = func_handles{2};
homds.Jacobian  = func_handles{3};
homds.JacobianP = func_handles{4};
homds.Hessians  = func_handles{5};
homds.HessiansP = func_handles{6};
homds.Der3 = func_handles{7};
homds.Der4 = func_handles{8};
homds.Der5 = func_handles{9};

cds.oldJac = [];
cds.oldJacX = [];
xp = [x;p(ap)];
cds.ndim = length(xp);

homds.x0 = x(1:homds.nphase);
  
pcell = num2cell(p);
%
% 3. COMPUTE A, A1
%--------------------------------------------
A = cjac(homds.func,homds.Jacobian,homds.x0,pcell,homds.ActiveParams);
A1 = cjacp(homds.func,homds.JacobianP,homds.x0,pcell,ap);
%
% 4. COMPUTE B, B1, B2, q0, q1, p0 and p1
%--------------------------------------------
% i)  B & B1
B = chess(homds.func,homds.Jacobian,homds.Hessians,homds.x0,pcell,homds.ActiveParams);
B1 = chessp(homds.func,homds.Jacobian,homds.HessiansP,homds.x0,pcell,homds.ActiveParams);
if size(B1,3) > length(ap)
    B1 = B1(:,:,ap);
end
%
% ii) B2
for i=homds.ActiveParams
    p1 = pcell; p1{i} = p1{i}-cds.options.Increment;
    p2 = pcell; p2{i} = p2{i}+cds.options.Increment;
    Hjp2 = cjacp(homds.func,homds.JacobianP,homds.x0,p2,ap); 
    Hjp1 = cjacp(homds.func,homds.JacobianP,homds.x0,p1,ap); 
    tmph(:,:,i) = Hjp2 - Hjp1;
end
tmph = tmph/(2*cds.options.Increment);
B2 = tmph;
if size(B2,3) > length(ap)
    B2 = B2(:,:,ap);
end
% iii) q0, q1, p0 & p1
  [X,D] = eig(A);
  index1 = find(abs(diag(D)) < 1e-3); %If ok, index1 is 1x2 array otherwise
  vext = real(X(:,index1(1)));
  [X,D] = eig(A');
  index1 = find(abs(diag(D)) < 1e-3);
  wext = real(X(:,index1(1)));
  Bord = [ A wext; vext' 0];
  bunit=[zeros(homds.nphase,1);1];
  q0=Bord\bunit; 
  q0=q0(1:homds.nphase);          % A q0 = 0 , <vext,q0> = 1
  p1=Bord'\bunit;
  p1=p1(1:homds.nphase);          % A'p1 = 0 , <wext,p1> = 1
  Bord = [ A p1; q0' 0];
  q1 = Bord\[q0; 0];		
  q1 = q1(1:homds.nphase);		  % A q1 = q0, <q0, q1>  = 0
  p0 = Bord'\[p1; 0];
  p0 = p0(1:homds.nphase);		  % A'p0 = p1, <p0, p1>  = 0
%
% 5. NORMALIZE SO THAT <p0,q0>=<p1,q1>=1 & <p0,q1>=<p1,q0>=0 
%--------------------------------------------
  mu = sqrt(abs(q0'*q0));
  q0 = (1/mu)*q0;
  q1 = (1/mu)*q1;
  q1 = q1 - (q0'*q1)*q0;
  nu = q0'*p0;
  p1 = (1/nu)*p1;
  p0 = p0 - (p0'*q1)*p1;
  p0 = (1/nu)*p0;
%
% 6. BT NORMAL FORM COEFFICIENTS (a,b)
%--------------------------------------------
hessIncrement = (cds.options.Increment)^(3.0/4.0);
if (cds.options.SymDerivative >= 2)
  hess = chess(homds.func,homds.Jacobian,homds.Hessians,homds.x0,pcell,homds.ActiveParams);
else
  hess = [];
end
M1 = multilinear2(homds.func,hess,q0,q0,homds.x0,pcell,hessIncrement);	% B(q0,q0)
M2 = multilinear2(homds.func,hess,q0,q1,homds.x0,pcell,hessIncrement);	% B(q0,q1)
M3 = multilinear2(homds.func,hess,q1,q0,homds.x0,pcell,hessIncrement);	% B(q1,q0)
M4 = multilinear2(homds.func,hess,q1,q1,homds.x0,pcell,hessIncrement);	% B(q1,q1)
% 
a = 0.5*p1'* M1;
b = p1'*M2+p0'*M1;
%
% 7. COMPUTE H200, H201, H01, K1, K2, H12 and H22 (THE HOMOCLINIC PREDICTOR COEFFICIENTS)
%--------------------------------------------
 for i=1:homds.nphase
    M5(:,i) = p1' * B(:,:,i) * q1;
    M6(:,i) = p1' * B(:,:,i) * q0;
    M7(:,i) = p0' * B(:,:,i) * q0;
    M8(:,i) = p0' * B(:,:,i) * q1;    
end
for i=1:length(ap)
    M9(:,i) = p1' * B1(:,:,i) * q1;
    M10(:,i) = p1' * B1(:,:,i) * q0;
    M11(:,i) = p0' * B1(:,:,i) * q0;
    M12(:,i) = p0' * B1(:,:,i) * q1;  
end
%
H200 = Bord\[2*a*q1-M1; 0];
H200 = H200(1:end-1,:);
H201 = Bord\[b*q1-M2+H200; 0];
H201 = H201(1:end-1,:);
%
C=[M6; M7+M5]            ; D=[M10; M11+M9];
F=[q1 zeros(size(A,1),1)]; G=[0.5*p1'*M4 0; -p0'*M4+3*p0'*H201 1];
tmpHK = [A A1; C D];
HK    = tmpHK\[F; G];
%
H01 = HK(1:end-2,:);
K1  = HK(end-1:end,:);
%
for i = 1:homds.nphase
    N1(:,i)  =  p1' *  B(:,:,i) * H01(:,2);
end
for i=1:length(ap)
    N2(:,i) =  p1' * B1(:,:,i) * H01(:,2);
    N3(:,i) =  p1' * B2(:,:,i) * K1(:,2) ;
end
N1  = N1 * H01(:,2);
N2 = N2 * K1(:,2);
N3 = N3 * K1(:,2);
K2   = -(N1 + 2 * N2 + N3) * K1(:,1);
%--------------------------------------------------------------------------
M13 = multilinear2(homds.func,hess,q0,H01(:,2),homds.x0,pcell,hessIncrement);	% B(q0,H01.1)
%
for i=1:length(ap)
    M14(:,i) =  q0' * B1(:,:,i); 
end
%
M14 = M14*K1(:,2);   % B1(q0,K1.1)
M15 = multilinear2(homds.func,hess,H01(:,2),H01(:,2),homds.x0,pcell,hessIncrement);	% B(H01.1,H01.1)
for i=1:length(ap)
    M16(:,i) =  H01(:,2)' * B1(:,:,i); 
    M17(:,i) =  B2(:,:,i) * K1(:,2);
end
M16 = M16*K1(:,2);   %B1(H01.1,K1.1)
M17 = M17*K1(:,2);   %B2(K1.1,K1.1)
z = M15 + 2*M16 + M17;
%
H12 = Bord\[-(M13+M14); 0];
H12 = H12(1:end-1,:);
%
H22 = Bord\[-(z+A1*K2); 0];
H22 = H22(1:end-1,:);
%--------------------------------------------------------------------------
%
% 8. INITIAL CYCLE
%--------------------------------------------
% a) THE INITIAL  APPROXIMATION OF ALPHA  
%
alpha = ((10*b)/(7*a)) * eps^2 * K1(:,2) + eps^4*((-4/a)*K1(:,1) + ...
    K1(:,2)*((b^3)/(a^3))*(288/2401) +((50*b^2)/(49*a^2))*K2);
p(ap) = p(ap) + alpha;  % tmpfreep = tmpfreep + alpha;
homds.P0 = p;
%
% b) THE INITIAL HALF-RETURN TIME VALUE T:-
%
homds.x0 = homds.x0 + eps^2 * ((10*b)/(7*a) * H01(:,2) + 2/a * q0);

x0 = homds.x0;
tmpmiddle = eps^2 * ((10*b)/(7*a) * H01(:,2) - 4/a * q0) + homds.x0;
%|w0(-inf,eps)-w0(T,eps)|=delta0 => (6*eps^2/|a|)*tanh(eps*T)^2=delta0
delta0 = norm(tmpmiddle - x0)/kfactor;
myarg = sqrt((delta0 * abs(a)) / (6 * eps^2));
myargn = log(1/myarg + ((1/myarg)^2 - 1));
if ~isempty(imag(myargn)) && (imag(myargn) ~= 0)
    myargn = log(1/myarg - ((1/myarg)^2 - 1));
end
homds.T = myargn/eps;
%
% c) THE INITIAL APPROXIMATION OF CYCLE
%
for i=1:length(homds.finemsh)    
    % CONVERSION FROM [0,1] --> [-T,+T]
    t = (2*homds.finemsh(i) - 1) * (homds.T);
    x0 = homds.x0;
    %
    ups(:,i) = eps^2 * ( (10*b)/(7*a) * H01(:,2) + (1/a) * (2 - 6 * sech(t * eps)^2) * q0) + ...
        eps^3*(1/a)*(12 * (sech(t * eps))^2 * tanh(t * eps) * q1 + (-(72/7)*sinh(t * eps)*b*...
        log(cosh(t * eps))/(cosh(t * eps)^3*a))*q0) + eps^4*((1/a)*(-72*b/(7*a)*(sinh(t * eps)^2 + ...
        (1-2*sinh(t * eps)^2)*log(cosh(t * eps)))/cosh(t * eps)^4)*q1 + ...
        (1/a)*(-(216/49)*log(cosh(t * eps))^2*(cosh(2*t*eps)-2)*b^2/(cosh(t * eps)^4*a^2)-...
        (216/49)*log(cosh(t * eps)) *(1 - cosh(2*t*eps))*b^2/(cosh(t * eps)^4*a^2)-...
        (18*b^2/49)*(6*t*eps *sinh(2*t*eps) - 7*cosh(2*t*eps) + 8)/(cosh(t * eps)^4*a^2))*q0 +...
        1/(2*a^2) * H200* ((2 - 6 * sech(t * eps)^2))^2 + H12 * 10*b/(7*a^2)*...
        ((2 - 6 * sech(t * eps)^2)) + (-4/a)*H01(:,1) + H01(:,2)*((b^3)/(a^3))*(288/2401) + H22*(50*b^2/(49*a^2)));
    %
    ups(:,i) = ups(:,i) + homds.x0;
end
%
%plot(ups(1,:),ups(2,:))
%
%d) THE INITIAL epsilon0 & epsilon1 :
%
homds.eps0 = norm(ups(:,1) - homds.x0);
homds.eps1 = norm(ups(:,end) - homds.x0);
%
% 9. DIMENSION OF STABLE & UNSTABLE SUBSPACE 
%--------------------------------------------
Hom_calc_weights;
A = cjac(homds.func,homds.Jacobian,x0,num2cell(p),homds.ActiveParams);
D = eig(A);
homds.nneg = sum(real(D) < 0);
homds.npos = homds.nphase-homds.nneg;
homds.Ysize = homds.nneg*homds.npos;
%
% 10. COMPOSE x
%--------------------------------------------
% 
% i. CYCLE
%
x1 = reshape(ups,size(ups,1)*size(ups,2),1);
v = []; 
[x1,v]=Hom_new_mesh(x1,v,ntst,ncol);
%
% ii. EQUILIBRIUM COORDINATES
%
x1 = [x1; x0];
%
% iii. TWO FREE PARAMETERS
%
x1 = [x1; homds.P0(homds.ActiveParams)];
%
% iv. EXTRA FREE PARAMETER (FREE HOMOCLINIC PARAMETER)
%
extravec = [homds.T; homds.eps0; homds.eps1];
x1 = [x1; extravec(find(homds.extravec))];

for i=1:homds.nneg
    x1 = [x1; zeros(homds.npos,1)];
end

for i=1:homds.npos
    x1 = [x1; zeros(homds.nneg,1)];
end
%
%
% 11. ASSIGN SOME VALUES TO HOMOCLINIC FIELDS
%--------------------------------------------
% a) YS AND YU, INITIALIZED TO 0
%
homds.YS = zeros(homds.npos,homds.nneg);
homds.YU = zeros(homds.nneg,homds.npos);

% b) THIRD PARAMETER = UNSTABLE_FLAG:-
%
% i) 1 IF WE WANT THE UNSTABLE SPACE. 
% ii) 0 IF WE WANT THE STABLE ONE.
[QU, se] = computeBase(A,0,homds.npos);
[QS, se] = computeBase(A,1,homds.nneg);

homds.oldStableQ = QS;
homds.oldUnstableQ = QU;
homds.ups = [];
homds.ndim = length(x1);
cd.ndim = homds.ndim;
x = x1;
%
%
% 12. INITIAL TANGENT VECTOR
%--------------------------------------------
ups0 = reshape(x(1:size(ups,1)*size(ups,2),:),homds.nphase,homds.tps);
pp1 = num2cell(p);
 for i=1:homds.tps
    homds.upoldp(:,i) = 2*T*feval(homds.func, 0, ups0(:,i), pp1{:});
 end
BTHomjac = BVP_Hom_jac(homds.func,x(1:size(ups,1)*size(ups,2),1),x0,p,homds.T,homds.eps0,homds.eps1,homds.YS,homds.YU);
[Q,R] = qr(full(BTHomjac)');
v = full(Q(:,end));

%---------------------------------------------------------------------------
function init_homds(odefile,x,p,ap,ntst,ncol,extravec,T,eps0,eps1,nph)
global homds 

homds.odefile = odefile;
func_handles = feval(homds.odefile);
homds.func = func_handles{2};
homds.Jacobian  = func_handles{3};
homds.JacobianP = func_handles{4};
homds.Hessians  = func_handles{5};
homds.HessiansP = func_handles{6};
homds.Der3=[];
siz = size(func_handles,2);

if siz > 9
    j=1;
    for k=10:siz
        homds.user{j}= func_handles{k};
        j=j+1;
    end
else homds.user=[];
end

homds.nphase = nph;
homds.ActiveParams = ap;
homds.P0 = p;
homds.extravec = extravec;

Hom_set_ntst_ncol(ntst,ncol,(0:ntst)/ntst);

homds.T = T;
homds.eps0 = eps0;
homds.eps1 = eps1;
homds.cols_p1 = 1:(homds.ncol+1);
homds.cols_p1_coords = 1:(homds.ncol+1)*homds.nphase;
homds.ncol_coord = homds.ncol*homds.nphase;
homds.col_coords = 1:homds.ncol*homds.nphase;
homds.pars = homds.ncoords+(1:3);
homds.phases = 1:homds.nphase;
homds.ntstcol = homds.ntst*homds.ncol;
homds.wp = kron(homds.wpvec',eye(homds.nphase));
homds.pwwt = kron(homds.wt',eye(homds.nphase));
homds.pwi = homds.wi(ones(1,homds.nphase),:);

homds.bialt_M1 = [];
homds.bialt_M2 = [];
homds.bialt_M3 = [];
homds.bialt_M4 = [];
homds.multipliers = nan;
homds.monodromy = [];
homds.multi_r1 = [];
homds.multi_r2 = [];
homds.ups = [];
homds.vps = [];
homds.tsts = 1:homds.ntst;
homds.cols = 1:homds.ncol;

homds.HTPstep = 0;