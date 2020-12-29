clear all
%%
% Degradation rate:
ku34 = 0.05;   kms = 0.5;   ks = 0.125;  ku200 = 0.05;   kmz = 0.5;   kz = 0.1;  ke=0.1;
kn=0.1;
% Transcription rate:
gu34 = 1350;   gms = 90;   gs = 100; gu200 = 2100;   gmz = 11;   gz = 100; ge=5000; gn=80000;
% Hills function threshold :
s0u34 = 300000;   s0ms = 200000;   z0u34 = 600000; u034 = 10000;  e0mz=20000; z0e=100000;
I0ms=50000; z0u200 = 220000;   z0mz = 25000;   s0u200 = 180000;   s0mz = 180000; u2000 = 10000;
n0z=800000; n0n=100000; n0e=100000; n0u200=500000;
% Cooperativity:
nsu34 = 1;    nsms = 1;   nsmz = 2;   nu34 = 2; nI = 2; nze=2; nnu200=4;
nzu200 = 3;   nsu200 = 2;   nzmz = 2;   nu200 = 6; nzu34=2; nemz=2; nnz=2; nnn=3; nne=5;
% fold change
lamdasu34 =0.1;   lamdazu34 = 0.2;  lamdasms = 0.1;   lamdaIms = 10; lamdaze =0.1;
lamdazu200 =0.1;   lamdasu200 = 0.1;  lamdazmz = 7.5;   lamdasmz = 10; lamdaemz=0.8;
lamdanz=2; lamdann=7; lamdane=8; lamdanu200=4;
% external signal
s = 329000;
u201=0;
u200=0;
mz2=1500;
N=1;

for i=1:20000;
%Ym function components
Mu0=1/(1+u200/u2000)^nu200;
Mu1=(u200/u2000)/(1+u200/u2000)^nu200;
Mu2=(u200/u2000)^2./(1+u200/u2000)^nu200;
Mu3=(u200/u2000)^3./(1+u200/u2000)^nu200;
Mu4=(u200/u2000)^4/(1+u200/u2000)^nu200;
Mu5=(u200/u2000)^5./(1+u200/u2000)^nu200;
Mu6=(u200/u2000)^6/(1+u200/u2000)^nu200;
%ZEB protein
z=(gz*mz2*(Mu0+0.6*6*Mu1+0.3*15*Mu2+0.1*20*Mu3+0.05*15*Mu4+0.05*6*Mu5+0.05*Mu6))/kz;

%Hills functions
Hillszu200=(1+lamdazu200*(z/z0u200)^nzu200)./(1+(z/z0u200)^nzu200);
Hillssu200=(1+lamdasu200*(s/s0u200)^nsu200)/(1+(s/s0u200)^nsu200);
Hillszmz=(1+lamdazmz*(z/z0mz)^nzmz)/(1+(z/z0mz).^nzmz);
Hillssmz=(1+lamdasmz*(s/s0mz)^nsmz)/(1+(s/s0mz)^nsmz);

Hillsze = (1+lamdaze*(z/z0e)^nze)/(1+(z/z0e)^nze);
Hillsnz = (1+lamdanz*(N/n0z)^nnz)/(1+(N/n0z)^nnz);
Hillsnn = (1+lamdann*(N/n0n)^nnn)/(1+(N/n0n)^nnn);
Hillsne = (1+lamdane*(N/n0e)^nne)/(1+(N/n0e)^nne);
Hillsnu200=(1+lamdanu200*(N/n0u200)^nnu200)/(1+(N/n0u200)^nnu200);
%ECAD
E=ge*Hillsze*Hillsne/ke;
Hillsemz = (1+lamdaemz*(E/e0mz)^nemz)/(1+(E/e0mz)^nemz);
%NFATcff
N=gn*Hillsnn/kn;
%Equations
mz1=(gu200*Hillszu200*Hillssu200*Hillsnu200-ku200*u200)/(0.005*6*Mu1+2*0.05*15*Mu2+3*0.5*20*Mu3+4*0.5*15*Mu4+5*0.5*6*Mu5+6*0.5*Mu6);
mz2=gmz*Hillszmz*Hillssmz*Hillsemz*Hillsnz/((0.04*6*Mu1+0.2*15*Mu2+20*Mu3+15*Mu4+6*Mu5+Mu6)+kmz);
amat(i,:)=[u200 mz1 mz2];
u200=u200+3;

end
plot(amat(:,1),amat(:,2),'b','LineWidth',2)
ylim([0,2500])
hold on
plot(amat(:,1),amat(:,3),'r','LineWidth',2)