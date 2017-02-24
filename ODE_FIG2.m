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
%ODE CODE OF FIGURE2
function y= ODE_FIG2(t,C)
global EA;
global EH;
global EF;
global ED;
global EG;
global E1;
global E2;
global E3;
global E4;
global E5;
global E6;
global K1;
global K2;
global K3;
global K4;
global K5;
global K6;
global K;


N = size(C);
%
for i=1:N
    if C(i) <= 0.0
			C(i) = 0;
	end
end
%}
Axt = C(1);
A = C(2);
B = C(3);
H = C(4);
D = C(5);
E = C(6);
F = C(7);
G = C(8);
VA = Axt*EA/(Axt+K);
V1 = A*E1/(A+K1);
V2 = B*E2/(B+K2);
V3 = B*E3/(B+K3);
V4 = B*E4/(B+K4);
V5 = E*E5/(E+K5);
V6 = E*E6/(E+K6);
VH = H*EH/(H+K);
VD = D*ED/(D+K);
VF = F*EF/(F+K);
VG = G*EG/(G+K);



y = [0;VA-V1;V1-V2-V3-V4;V2-VH;V3-VD;V4-V5-V6;V5-VF;V6-VG];

