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

%HCSA CODE OF FIGURE3
%IN THIS CODE:
%"LP" MEANS Linear Programming
%"ECS" IS AN ID, DOESN'T HAVE ANY MEANINGS
%"_eq"，"_ieq" MEANS INFORMATION IN THIS MATRIX WILL BE USED AS A EQUITY/INEQUITY CONSTRANT IN HCSA

clear
clc
figure

%%数据点的编号

thrsh = 0.1;
MaxNum = 10;

load('FIG3_ECS_eq_10.mat');
ECS_eq_ori = ECS_eq;

load('FIG3_ECS_ieq_10.mat');
load('FIG3_E_vary_whole.mat');
load('FIG3_norm.mat')
load('FIG3_experimental.mat')

%INITIALIZATION OF LINEAR PROGRAMMING
M = 100;
beq = zeros(4,1);
f = [zeros(8,1);-1];
b1 = [zeros(6,1);ones(6,1)];
b2 = [zeros(3,1);ones(3,1)];
lb=zeros(9,1);
ub=M*ones(9,1);

%DATUM EV USED IN FIG3
Iori = [100,100,100,60 48,1584.89319246111,28,6];
                         
%191 TRANSCRIPTIONAL LEVEL CHANGE, HCSA GIVES 191 PREDICTIONS
for a =1:191
	beq = zeros(4,1);
	ECS_eq = ECS_eq_ori;
	set =a;
	I= Iori.*E_vary(set,:);
	%I
	I = I/max(I);
	
    %WE CUT THE REACTIONS NETWORK(FIG3) INTO 2 NETWORKS
	%CENTERED METABOLITE PDVA AND PVA, RESPECTIVELY
	
	%DOING HCSA ONTO THE ONE CENTERED PDVA
	I1 = I(1:6);
	I1 = I1/max(I1);
	ECS_1_ieq(1:6,9) = I1';
    ECS_1_ieq(7:end,9) = ones(6,1)-I1';
    ECS_see = ECS_eq;
	%HCSA
    x0 = linprog(f,ECS_1_ieq,b1,ECS_eq,beq,lb,ub);
    l =x0;
    x0 = x0/norm(x0(1:6));
                        
    %INITIALIZATION OF FLAG
	flag = 0;
	count = 0;
	prin = I1./norm(I1);
	while flag == 0
		error = 0;
		if count>MaxNum
			break;
		end
		for i =1:6
			%IF EV AND FV DON'T CORELATED WELL, WE CHANGE THE UPPER BOUND OF FLUXES------------------
			%THIS MODIFICATION OF UPPER BOUND MAY LEAD THE HYPER-CUBE SHRINK INTO A SMALLER VOLUNM---
			%A BETTER OPTIMIZATION METHOD COULD REPLACE MY NAIVE ONE, IF YOU HAVE ANY IDEA, CONTACT ME----
			if abs(log2(x0(i)./prin(i)'))>thrsh
				ECS_1_ieq(6+i,9) = max(ECS_1_ieq(6+i,9)+(1/(count+1))*log2(x0(i)/prin(i)),0);
				error = 1;
			end
			%----------------------------------------------------------------------------------------
		end
		x0 = linprog(f,ECS_1_ieq,b1,ECS_eq,beq,lb,ub);
		l = x0;
		x0 = x0/norm(x0(1:6));
		count = count+1;
		if error == 0
			flag = 1;
		end            
	end
    lpr1 = x0;
    lpr1save(a,:) = lpr1;
    ECS_eq(5,:) = [lpr1(2) -lpr1(1) 0 0 0 0 0 0 0];
    ECS_eq(6,:) = [lpr1(3) 0 -lpr1(1) 0 0 0 0 0 0];
    ECS_eq(7,:) = [lpr1(4) 0 0 -lpr1(1) 0 0 0 0 0];
    ECS_eq(8,:) = [lpr1(5) 0 0 0 -lpr1(1) 0 0 0 0];
    ECS_eq(9,:) = [lpr1(6) 0 0 0 0 -lpr1(1) 0 0 0];
    beq = zeros(9,1);
    
	%DOING HCSA ONTO THE ONE CENTERED PVA
	I2 = I(6:8);
    I2 = I2/max(I2);
    ECS_2_ieq(1:3,9) = I2';
    ECS_2_ieq(4:end,9) = ones(3,1)-I2';
    x0 = linprog(f,ECS_2_ieq,b2,ECS_eq,beq,lb,ub);
    x0 = x0/norm(x0(6:8));
	
	%INITIALIZATION OF FLAG	
	flag = 0;
	count = 0;
	prin = I2./norm(I2);
	while flag == 0
		error = 0;
		if count>MaxNum
			break;
		end
		for i =1:3
			%IF EV AND FV DON'T CORELATED WELL, WE CHANGE THE UPPER BOUND OF FLUXES------------------
			%THIS MODIFICATION OF UPPER BOUND MAY LEAD THE HYPER-CUBE SHRINK INTO A SMALLER VOLUNM---
			%A BETTER OPTIMIZATION METHOD COULD REPLACE MY NAIVE ONE, IF YOU HAVE ANY IDEA, CONTACT ME----
			if abs(log2(x0(5+i)./prin(i)'))>thrsh
				ECS_2_ieq(3+i,9) = max(ECS_2_ieq(3+i,9)+(1/(count+1))*log2(x0(5+i)/prin(i)),0);
				error = 1;
			end
		end
		x0 = linprog(f,ECS_2_ieq,b2,ECS_eq,beq,lb,ub);
		x0 = x0/norm(x0(6:8));
		count = count+1;
		if error == 0
			flag = 1;
		end            
	end
	x0 = linprog(f,ECS_2_ieq,b2,ECS_eq,beq,lb,ub);
    lpr2 = x0;
    lpr = [lpr2(4) lpr2(5) lpr2(7) lpr2(8)];
    lpr = lpr/norm(lpr);
    lprsave(a,:) =  lpr;
end

%pred_test_whole RECORDS THE NORM(ABSOLUTE VALUE) OF FV
%OBTAINED BY MULTIPLE LINEAR REGRESSION:
% NORM_OF_FV = f(DOSE OF ENZYME)
%TRAINING SET INDEX: 1 3 6 9 15 17 18 19 20 21 22 23 24 25;
for i = 1:191
	lprsave(i,:) = lprsave(i,:)./norm(lprsave(i,:));
	lprsave(i,:) = lprsave(i,:).*pred_test_whole(i);
end
%PLOT FIGURE3(B,C,D,E)
loglog(lprsave,Exp,'.');
