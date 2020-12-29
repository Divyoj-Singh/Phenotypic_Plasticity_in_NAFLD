clear
global x v s h f opt
init

OPTIONS = [];
hls = MyML;
[t,y] = ode45(hls{2},[0 1000],[0 0],OPTIONS,30,10);
x0 = y(end,:)';

opt=contset;opt=contset(opt,'Singularities',1);
opt = contset(opt,'MaxStepsize',10);
opt=contset(opt,'MaxNumpoints',2500);
[x1,v1]=init_EP_EP(@MyML,x0,[30;6],[1]);
disp('>> [x,v,s,h,f]=cont(@equilibrium,x0,[],opt);');
[x,v,s,h,f]=cont(@equilibrium,x1,[],opt);

x1=x(1:2,s(5).index);p=[x(end,s(5).index);6];
[x0,v0]=init_H_LC(@MyML,x1,p,[1],1e-6,40,4);
opt = contset(opt,'MaxNumPoints',200);
opt = contset(opt,'Multipliers',0);
opt = contset(opt,'Adapt',1);
opt = contset(opt,'MaxStepsize',5);
opt = contset(opt,'FunTolerance',1e-6);
opt = contset(opt,'VarTolerance',1e-6);
disp('>> [xlc,vlc,slc,hlc,flc]=cont(@limitcycle,x0,v0,opt);');
[xlc,vlc,slc,hlc,flc]=cont(@limitcycle,x0,v0,opt);

[x0,v0]=init_LC_LC(@MyML,xlc,vlc,slc(end),[1],40,4);
opt = contset(opt,'PRC',1);
opt = contset(opt,'dPRC',1);
opt = contset(opt,'Input',1);
opt = contset(opt,'MaxNumPoints',10);
[xlc2,vlc2,slc2,hlc2,flc2]=cont(@limitcycle,x0,v0,opt);
