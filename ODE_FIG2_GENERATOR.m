%-------------------
%AUTHOR: TIANYU ZHANG
%EMAIL: zhangtianyu@pku.edu.cn
%PUBLICATION:
%	Z.W.,XIE, T.Y.,ZHANG, Q.OUYANG, Genome-scale fluxes were predicted 
%	under the guide of enzyme abundance using a novel 
%	Hyper-Cube Shrink Algorithm
%IF YOU'VE GOT ANY QUESTIONS ABOUT THIS .m FILE,
%FEEL FREE TO CONTACT ME.
%:)
%--------------------
clear
clear global all

%INITIALIZATION OF THE ENZYMATIC PARAMETERS
global EA;EA = 4;
global EH;EH = 1;
global EF;EF = 1;
global ED;ED = 1;
global EG;EG = 1;
global E1;E1ori = 4;
global E2;E2ori = 1;
global E3;E3ori = 1;
global E4;E4ori = 2;
global E5;E5ori = 1;
global E6;E6ori = 1;
global K;K = 1;
global K1;K1 = 1;
global K2; K2 = .25;
global K3; K3 = .5;
global K4; K4 = .7;
global K5; K5 = .7;
global K6; K6 = 1;
E1 = E1ori;
E2 = E2ori;
E3 = E3ori;
E4 = E4ori;
E5 = E5ori;
E6 = E6ori;

init = [1 0 0 0 0 0 0 0];
time = [0:500];
raaw = 1;

%CHANGE ONE SINGLE PARAMETER
%IN THIS EXAMPLE:E6
for i = -1:.02:1
    k = 10^(i);
    E6 = E6ori*k;
    [t,C]=ode15s(@ODE_FIG2, time, init);
	%C IS CONTAINS THE CONCENTRATION OF METABOLITES
Axt = C(end,1);
A = C(end,2);
B = C(end,3);
H = C(end,4);
D = C(end,5);
E = C(end,6);
F = C(end,7);
G = C(end,8);
%CALCULATE THE FLUX
V1 = A*E1/(A+K1);
V2 = B*E2/(B+K2);
V3 = B*E3/(B+K3);
V4 = B*E4/(B+K4);
V5 = E*E5/(E+K5);
V6 = E*E6/(E+K6);
V = [V1 V2 V3 V4 V5 V6];
%ODE6 CONTAINS THE FLUX DISTRIBUTION EDITING ENZYME 6
ODE6(raaw,1:6)  = V;
raaw = raaw+1;
end
%PLOT FIG2(I)
xaxis = -1:.02:1;
xaxis = 10.^(xaxis);
figure;plot(xaxis,ODE6);
