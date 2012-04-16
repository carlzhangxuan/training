%matlab code-[dynare tool box]-solow T case
%----------------------------------------------------------------
%G->0.1->1.1G  by T/Sita = 0.2262 /gamma =0.5475
%z->1.1z
%----------------------------------------------------------------
close all;
%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------
var w Y K C_T G A L I KG q Pi m R F Z Tl;
varexo z;
parameters a b Sita Sigma gamma T r zB zG phoL phoK ItaL ItaK delta EY mu TC mI;
%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------
b  = 0.018;
a = 0.0388;
delta = 0.1;
EY = 0.29;
r = 0.04;
zB = 0.532;
zG = 0.532;
Sigma = 0.5;
phoL = 1.1;
phoK = 0.8;
ItaL = 0;
%(Solow-neutral case:0;Harrod-neutral case:0.12;Hicks-neutral case:0.08)-SigmaY=0.5
ItaK = 0.29;
%(Solow-neutral case:0.29;Harrod-neutral case:0;Hicks-neutral case:0.08)-SigmaY=0.5
%----------handle------------------------------------------------
Sita = 0.2262; 
gamma =0.5475;
TC = 0.01;
mu = 0.01;
%Tl =1;
mI =0.001;
T = 0.01 
%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------
model;
G = mu*m+Tl*T+TC*C_T;
G = z*0.0126826;
%----------handle------------------------------------------------
K-K(-1)=(zB*ln(((I/K(-1))+zB)/zB)-delta)*K(-1);
KG-KG(-1)=(zG*ln(((G/KG(-1))+zG)/zG)-delta)*KG(-1);
q(+1)-q =q*(r-zB*ln(((I/K(-1))+zB)/zB)+delta)+(I/K)-EY*(phoK*KG^(ItaK))^((Sigma-1)/Sigma)*(Y/K)^(1/Sigma); 
A-A(-1)=r*A(-1)+w*L-C_T*(1+TC)-Tl*T-(Pi+r)*m;
%TC(+1)->TC
C_T(+1)-C_T=(r-a-((TC-TC)/(1+TC)))*C_T-b*Sita*(a+b)*(A+(Pi+r)*m)/((gamma+Sita)*(1+TC));
Pi = R-r;
m-m(-1)=(mu-Pi)*m;
1 = q*(zB/(zB+(I/K(-1))));
Y = (EY*(phoK*KG^(ItaK)*K)^((Sigma-1)/Sigma)+(1-EY)*(phoL*KG^(ItaL)*L)^((Sigma-1)/Sigma))^(Sigma/(Sigma-1));
w = (1-EY)*(phoL*KG^ItaL)^((Sigma-1)/Sigma)*(Y/L)^(1/Sigma);
(1-L)=gamma*C_T*(1+TC)/(w*Sita);
%m(+1) =(1-Sita-gamma)*C_T*(1+TC)/(R*Sita)-b*m;
m-m(-1)=(r-a-((R(+1)-R)/R))*m-b*m+mI;
A=q*K+m+F;
Z=Y-I-C_T-G;
end;
%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------
initval;
  G = 1;
  C_T = 0.436787;
  K = 4.75617;
  Y = 2.25692;
  Pi = 0.01;
  R = 0.0454788;
  m = 2;
  w = 3.12765;
  L = 0.480896;
  z = 1;  
  KG = 2.99956;
end;
steady;
check;
endval;
z = 1.1;
end;
steady;
check;
shocks;
var z;
periods 0:0;
values 1.1;
end;
simul(periods = 200);
%----------------------------------------------------------------
% 5. Some Results
%----------------------------------------------------------------

matlab code-[data]
clear all;
dynare SolowfixedG_T;
save Y_T.dat Y -ascii;
save q_T.dat q -ascii;
save K_T.dat K -ascii;
save L_T.dat L -ascii;
save C_T.dat C_T -ascii;
save Tl_T.dat Tl -ascii;

matlab code-[output]
clear all;

load Y_T.dat;load Y_Tc.dat;load Y_Mu.dat;
load q_T.dat;load q_Tc.dat;load q_Mu.dat;
load K_T.dat;load K_Tc.dat;load K_Mu.dat;
load L_T.dat;load L_Tc.dat;load L_Mu.dat;
load C_T.dat;load C_Tc.dat;load C_Mu.dat;
load F_T.dat;load F_Tc.dat;load F_Mu.dat;

figure;
subplot(321);plot((Y_T-Y_T(1))/Y_T(1),'-');holdon;plot((Y_Tc-Y_Tc(1))/Y_Tc(1),'--');hold on;plot((Y_Mu-Y_Mu(1))/Y_Mu(1),'-.');
axis([2 200 -0.01 0.02]);
title('\fontname{arial}Y');xlabel('\fontname{arial}periods','fontsize',8);legend('T','TC','Mu');ylabel('\fontname{arial}Output','fontsize',8);set(gca,'fontname','arial');
subplot(322);plot((q_T-q_T(1))/q_T(1),'-');holdon;plot((q_Tc-q_Tc(1))/q_Tc(1),'--');hold on;plot((q_Mu-q_Mu(1))/q_Mu(1),'-.');axis([2 200 -0.005 0.001]);
title('\fontname{arial}q');xlabel('\fontname{arial}periods','fontsize',8);ylabel('\fontname{arial}q','fontsize',8);set(gca,'fontname','arial');
subplot(323);plot((K_T-K_T(1))/K_T(1),'-');holdon;plot((K_Tc-K_Tc(1))/K_Tc(1),'--');hold on;plot((K_Mu-K_Mu(1))/K_Mu(1),'-.');axis([1 200 -0.008 0.006]);
title('\fontname{arial}K');xlabel('\fontname{arial}periods','fontsize',8);ylabel('\fontname{arial}Capital Stock','fontsize',8);set(gca,'fontname','arial');
subplot(324);plot((L_T-L_T(1))/L_T(1),'-');holdon;plot((L_Tc-L_Tc(1))/L_Tc(1),'--');hold on;plot((L_Mu-L_Mu(1))/L_Mu(1),'-.');axis([2 200 -0.01 0.015]);
title('\fontname{arial}L');xlabel('\fontname{arial}periods','fontsize',8);ylabel('\fontname{arial}Labor','fontsize',8);set(gca,'fontname','arial');
subplot(325);plot((C_T-C_T(1))/C_T(1),'-');holdon;plot((C_Tc-C_Tc(1))/C_Tc(1),'--');hold on;plot((C_Mu-C_Mu(1))/C_Mu(1),'-.');axis([2 200 0.002 0.012]);
title('\fontname{arial}C');xlabel('\fontname{arial}periods','fontsize',8);ylabel('\fontname{arial}Consumption','fontsize',8);set(gca,'fontname','arial');
subplot(326);plot((F_T-F_T(1))/F_T(1),'-');holdon;plot((F_Tc-F_Tc(1))/F_Tc(1),'--');hold on;plot((F_Mu-F_Mu(1))/F_Mu(1),'-.');axis([1 200 -0.6 0.4]);
title('\fontname{arial}F');xlabel('\fontname{arial}periods','fontsize',8);ylabel('\fontname{arial}F','fontsize',8);set(gca,'fontname','arial');