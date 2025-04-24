

clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l1=0.074;
l2=0.074; 
l3=0.074;


% Input parameters 
k=7;
l=13;
m=25;
n=19;
q=1624.88;
T=5.63;

%%%%%%%%%%%%%%%%%%%%%%%%%

alpha_i=0;

beta_i=0;

tr_ki=0;

tp_ki=0;

u_i=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=50000;
tB=300;
lambda=160000;
sigma=140000;
alpha=0.02;
E_delta=0.75;
x=8500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%remanufacturing quantity in remanufacturing time
rbar_ik=(sigma/(1+l1))*((1+tr_ki+k*T)^(1+l1)-(1+tr_ki)^(1+l1));


%production quantity would have been produced if the disruption had not occurd 
rbar_im=(sigma/(1+l1))*((1+tr_ki+m*T)^(1+l1)-(1+tr_ki)^(1+l1));

%theoritrical time requierd to produce 
t_alpha_r_i=(((1+l1)/sigma)*(rbar_ik+alpha_i)+1)^(1/(1+l1))-1; 

%ratio
Cr_i=tB/t_alpha_r_i;

%forgetting function 
fr_i=(l1*(1-l1)*log(rbar_ik+alpha_i))/(log(Cr_i+1));


%remembered quantity 
alpha_i1=((rbar_ik+alpha_i)^((l1+fr_i)/l1))*((rbar_im+alpha_i)^((-fr_i)/l1));


%equivalent time
tr_ki=(((1+l1)/sigma)*alpha_i1+1)^(1/(1+l1))-1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%manufacturing quantity in manufacturing time
qbar_in=(lambda/(1+l2))*((1+tp_ki+(n-k)*T)^(1+l2)-(1+tp_ki)^(1+l2));


%production quantity would have been produced if the disruption had not occurd 
qbar_im=(lambda/(1+l2))*((1+tp_ki+m*T)^(1+l2)-(1+tp_ki)^(1+l2));

%theoritrical time requierd to produce 

t_beta_p_i=(((1+l2)/lambda)*(qbar_in+beta_i)+1)^(1/(1+l2))-1; 

%ratio
Cp_i=tB/t_beta_p_i;

%forgetting function 
fp_i=(l2*(1-l2)*log(qbar_in+beta_i))/(log(Cp_i+1));


%remembered quantity 
beta_i1=((qbar_in+beta_i)^((l2+fp_i)/l2))*((qbar_im+beta_i)^((-fp_i)/l2));


%equivalent time
tp_ki=(((1+l2)/lambda)*beta_i1+1)^(1/(1+l2))-1;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%inspection

%%%%expected number of collected units
E_L=alpha*(l*q+(m-l)*(1-E_delta)*q)+(m-l)*E_delta*q;

%E_L=40000;

%%% if interruption not occurd
V_i1=((1-l3)*m*T*x+(u_i)^(1-l3))^(1/(1-l3))-u_i;

%%inspection time for i-th cycle
tauc_i=(1/((1-l3)*x))*((u_i+E_L)^(1-l3)-u_i^(1-l3));

c_i=tB/tauc_i;

%%forgetting exponent
fc_i=(l3*(1-l3)*log(u_i+E_L))/log(c_i+1);

%%initial experience for next ccycle
u_i1=((u_i+E_L)^((l3+fc_i)/l3))*(u_i+V_i1)^(-fc_i/l3);


%%%%%%%%%%%%%%%%print%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('forfetitting rate f_ri = %f\n',fr_i)
fprintf('forfetitting rate f_pi = %f\n',fp_i)
fprintf('forfetitting rate f_ci = %f\n',fc_i)
fprintf('equivalent time tr_ki= %f\n',tr_ki)
fprintf('equivalent time for tp_ki= %f\n',tp_ki)
fprintf('equivalent quantity for inspection u_i1= %f\n',u_i1)
fprintf('total return quantity E_L= %f\n', E_L)
fprintf('The value of alpha_i1 = %f\n',alpha_i1)
fprintf('The value of beta_i1 = %f\n',beta_i1)
