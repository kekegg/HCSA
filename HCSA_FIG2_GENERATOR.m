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

%HCSA CODE OF FIGURE2
%IN THIS CODE:
%"LP" MEANS Linear Programming
%"ECS" IS AN ID, DOESN'T HAVE ANY MEANINGS
%"_eq"£¬"_ieq" MEANS INFORMATION IN THIS MATRIX WILL BE USED AS A EQUITY/INEQUITY CONSTRANT IN HCSA

clear
clc

%INITIALIZATION OF LINEAR PROGRAMMING
load('FIG2_ECS_eq.mat');
ECS_eq_ori = ODE_ECS_eq;
load('FIG2_ECS_2_ieq.mat');
load('FIG2_ECS_1_ieq.mat');
M = 100;
f = [zeros(6,1);-1];
beq = zeros(2,1);
b2 = [zeros(3,1);ones(3,1)];
b1 = [zeros(4,1);ones(4,1)];
lb=zeros(7,1);
ub=M*ones(7,1);
%DATUM EV USED IN FIG2
Iori = [4,1.33333333333333,0.999999999999999,1.66666666666667,0.691573613268316,0.975093053398352];
raaw = 1;
thrsh = 0.1;
MaxNum = 10;
%CHANGE TRANSCRIPTIONAL LEVEL AND PREDICT FLUX DISTRIBUTION
for mul = -1:.02:1
    k = 10^(mul);
    beq = zeros(2,1);
    ECS_eq = ECS_eq_ori;
    I = Iori;
	%IN THIS EXAMPLE, EXPRESSION OF ENZYME 6 IS MODIFIED
    I(6) = I(6)*k;
    I = I/max(I);
	
	%WE CUT THE 6 REACTIONS NETWORK(FIG3) INTO 2 NETWORKS
	%CENTERED METABOLITE A AND B, RESPECTIVELY
	
	%DOING HCSA ONTO THE ONE CENTERED A
    I1 = I(1:4);
    I1 = I1/max(I1);
	%LOWER BOUND
    XIE_new_ECS_ieq_1(1:4,7) = I1';
	%UPPER BOUND
    XIE_new_ECS_ieq_1(5:end,7) = ones(4,1)-I1';
	%HCSA
    x0 = linprog(f,XIE_new_ECS_ieq_1,b1,ECS_eq,beq,lb,ub);
    x0 = x0/norm(x0(1:4));
    
	%INITIALIZATION OF FLAG
	flag = 0;
	count = 0;
	prin = I1;
	while flag == 0
		error = 0;
		if count>MaxNum
			break;
		end
		
		for i =1:4
			%IF EV AND FV DON'T CORELATED WELL, WE CHANGE THE UPPER BOUND OF FLUXES------------------
			%THIS MODIFICATION OF UPPER BOUND MAY LEAD THE HYPER-CUBE SHRINK INTO A SMALLER VOLUNM---
			%A BETTER OPTIMIZATION METHOD COULD REPLACE MY NAIVE ONE, IF YOU HAVE ANY IDEA, CONTACT ME----
			if abs(log2(x0(i)./prin(i)'))>thrsh
				XIE_new_ECS_ieq_1(4+i,7) = max(XIE_new_ECS_ieq_1(4+i,7)+(1/(count+1))*log2(x0(i)/prin(i)),0);
				error = 1;
			end
			%----------------------------------------------------------------------------------------
		end
		
		%HCSA
		x0 = linprog(f,XIE_new_ECS_ieq_1,b1,ECS_eq,beq,lb,ub);
		l = x0;
		x0 = x0/norm(x0(1:4));
		count = count+1;
		
		if error == 0
			flag = 1;
		end            
	end
	%%---------------------------------------------------------------------
	
	lpr1 = x0;
	lpr1save(raaw,:) = lpr1;
	ECS_eq(3,:) = [lpr1(2) -lpr1(1) 0 0 0 0 0 ];
	ECS_eq(4,:) = [lpr1(3) 0 -lpr1(1) 0 0 0 0];
    ECS_eq(5,:) = [lpr1(4) 0 0 -lpr1(1) 0 0 0 ];
    beq = zeros(5,1);
	
	%DOING HCSA ONTO THE ONE CENTERED B
	I2 = I(4:6);
	I2 = I2/max(I2);
	XIE_new_ECS_ieq_2(1:3,7) = I2';
	XIE_new_ECS_ieq_2(4:end,7) = ones(3,1)-I2';
	x0 = linprog(f,XIE_new_ECS_ieq_2,b2,ECS_eq,beq,lb,ub);
	x0 = x0/norm(x0(4:6));
							
	%INITIALIZATION OF FLAG				
	flag = 0;
	count = 0;
	prin = I2;
	while flag == 0
		error = 0;
		if count>MaxNum
			break;
		end
		for i =1:3
			%IF EV AND FV DON'T CORELATED WELL, WE CHANGE THE UPPER BOUND OF FLUXES------------------
			%THIS MODIFICATION OF UPPER BOUND MAY LEAD THE HYPER-CUBE SHRINK INTO A SMALLER VOLUNM---
			%A BETTER OPTIMIZATION METHOD COULD REPLACE MY NAIVE ONE, IF YOU HAVE ANY IDEA, CONTACT ME----
			if abs(log2(x0(3+i)./prin(i)'))>thrsh
				XIE_new_ECS_ieq_2(3+i,7) = max(XIE_new_ECS_ieq_2(3+i,7)+(1/(count+1))*log2(x0(3+i)/prin(i)),0);
				error = 1;
			end
		 end
		x0 = linprog(f,XIE_new_ECS_ieq_2,b2,ECS_eq,beq,lb,ub);
		x0 = x0/norm(x0(4:6));
		count = count+1;
		if error == 0
			flag = 1;
		end            
	end

    x0 = linprog(f,XIE_new_ECS_ieq_2,b2,ECS_eq,beq,lb,ub);
    lpr2 = x0;
    lpr = lpr2;
    noorm(raaw) = norm(lpr);
    lpr = lpr/norm(lpr);
    %LP6 CONTAINS THE FLUX DISTRIBUTION EDITING ENZYME 6
    LP6(raaw,:) =  lpr;
    raaw=raaw+1;
end

LP6 = LP6(:,1:6);
xaxis = -1:.02:1;
xaxis = 10.^(xaxis);
%PLOT FIG2(F)
figure;plot(xaxis,LP6);