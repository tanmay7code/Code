clc; 
clear; 
close all;

nv=6; % No of variables
lmt=[0,100000;0,100000;0,100000; 0, 100000; 0,10000; 0, 100000]; %lower and upper bound of first and second variables
d1=[2,2,2,2,2,2]; %accuracy of variables (number of digits after decimal place)
pop=100; %populations size
maxgen=100; %maximum no of generations to be performed
pc=0.8; %cross over probability
pm=0.01; %mutation probability

% Model parameters
d=50000; % demand rate 
s_m=400;  %manufacturing setup cost for the manufacturer for i-th cycle
d_m=120; %material cost for the manufacturer
gamma_m=24000; %manufacturing cost for the manufacturer
gamma_r=30000; %remanufacturing cost for the manufacturer
c_m=150;  %Warrenty cost for each defective inventory product by manufacturer
s_mr=25; %shipment transportation cost for the manufacture
sigma=140000; %initial rmanufacturing rate for the first cycle
lambda=160000; %initial manufacturing rate for the first cycle
hm1=3;  %remanufacturing inventory holding cost for the manufacturer
hm2=4;  %manufacturing inventory holding cost for the manufacturer
A_r=300; %order cost for the retailer
w_r=150; %purchasing cost for the retaile%
s_r_dash=15000; %inspection cost
I_r=175200; %inspection rate
h_r1=2; %retailer holding cost
h_r2=3; % retailer holding cost
h_c=5; % collction center holding cost
alpha=0.5; %percentage of used products returned by the custome
c_r=5; %offer price for the customer to return the product
beta=0.8; %percentage of the collected inventory that is collected from the customer and retailer at the collection center is eligible for remanufacturing;
x1=8500; %initial inspection rate of the collected inventory at the collection center

a_re=3*10^(-7);
b_re=0.0012;
c_re=1.4;
a_me=3*10^(-7);
b_me=0.0012;
c_me=1.4;
L_mr=100;
L_rc=100;
L_cm=100;

E=36.75; %vehicle fuel consumption (litre/km)
v_e=2.6; %vehicle emission standard
w_e=500; %warehouse emission standard for electricity
d_e=100; % emission for dfective
c_t=6;  %carbon tax
C_max=600000; % carbon limit
u=6;% carbon buying and selling price
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For the first cycle
E_L=0;

trk1=0;

tpk1=0;

u_1=0;

l_1=0.074;

l_2=0.074;

l_3=0.074;

E_delta=0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%binary_individual_form
% [slen,tlen]=binary_individual_form(nv,lmt,d);
slen=zeros(nv,1);
tlen=0; %Total length of an individual
for j=1:nv 
    ltmp =log((lmt(j,2)-lmt(j,1))*10^d1(j)+1)/log(2);%
    slen(j,1)=ceil(ltmp);%Total number elements for the jth variable
    tlen=tlen+slen(j,1);
end

% %binary_individual_initialize
% %binary_individual_decode
bp=zeros(pop,tlen);
xp=zeros(pop,nv);

for p=1:pop

    while (xp(p,1)>=xp(p,2)) || (xp(p,2)>=xp(p,4)) || (xp(p,4)>=xp(p,3))
        
        for s=1:tlen

        if rand<=0.5
            bp(p,s)=0;
        else
            bp(p,s)=1;
        end
         end

          c=0;
    for j=1:nv
        sum=0;
        for k=1:slen(j,1)
            sum=sum+bp(p,c+k)*(2^(slen(j,1)-k));
        end
        xp(p,j)=lmt(j,1)+sum*((lmt(j,2)-lmt(j,1))/(2^slen(j,1)-1));
        c=c+slen(j,1);
    end

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function_evalution 1
%function_evalution 2
% fp=user_funct(pop,xp);
%function value for parents

fp1=zeros(pop,1);
fp2=zeros(pop,1);

for p=1:pop

y1=xp(p,1);
y2=xp(p,2);
y3=xp(p,3);
y4=xp(p,4);
y5=xp(p,5);
y6=xp(p,6);

A1=0;
S1=zeros;
for j=1:ceil(y1)
    S1(j)=((1+trk1+j*y6)^(2+l_1)-(1+trk1+(j-1)*y6)^(2+l_1));
A1=A1+S1(j);
end

A2=0;
S2=zeros;
for j=1:ceil(y1)
    S2(j)=((1+trk1+(j-1)*y6)^(1+l_1));
A2=A2+S2(j);
end



A3=0;
S3=zeros;
for j=1:ceil(y1)-1
    S3(j)=((1+trk1+j*y6)^(1+l_1)-(1+trk1)^(1+l_1));
A3=A3+S3(j);
end

A4=0;
S4=zeros;
for j=ceil(y2)+1:ceil(y4)-1
    S4(j)=((1+tpk1+(j-y1+1)*y6)^(2+l_2)-(1+tpk1+(j-y1)*y6)^(2+l_2));
A4=A4+S4(j);
end



f1=@(y1,y2,y3,y4,y5,y6) (((d/(y2*y5+(y3-y2)*(1-E_delta)*y5))*(s_m+(y3-y2)*y5*d_m+(y3-y2)*y6*gamma_m+y2*y6*gamma_r+(y3-y2)*E_delta*y5*c_m+y3*y5*s_mr+(sigma/((1+l_1)*(2+l_2))*A1-(sigma*y6)/(1+l_1)*A2+(sigma*y6)/(1+l_1)*A3+(y2-y1)*(sigma*y6)/(1+l_1)*((1+trk1+y1*y6)^(1+l_1)-(1+trk1)^(1+l_1))-1/2*y2*(y2-1)*y5*y6)*hm1+(lambda/((1+l_2)*(2+l_2))*((1+tpk1+(y2-y1+1)*y6)^(2+l_2)-(1+tpk1)^(2+l_2))-1/(1+l_2)*(y2-y1+1)*lambda*y6*(1+tpk1)^(1+l_2)+lambda/((1+l_2)*(2+l_2))*A4-1/2*(y4^2-y4-y2^2+y2)*(1+tpk1)^(1+l_2)+y6*(y3-y4)*lambda/(1+l_2)*((1+tpk1+(y4-y1)*y6)^(1+l_2)-(1+tpk1)^(1+l_2))-1/2*(y3-y2-1)*y5*y6)*hm2+A_r+y3*y5*w_r+(y3-y2)*y5*s_r_dash*1/I_r+1/2*y3*(y3+1)*E_delta*y5*y6*h_r1+(1/2*y2^2*(y5-d*y6)*y6+1/2*y2*y5*y6+y2*(y3-y2)*y6*(y5-d*y6)+1/2*(y3-y2+1)*(y3-y2)*y6*((1-E_delta)*y5-d*y6)+(y3-y2)*y5*y6-(y3-y2)*(y6-y5/I_r)*(E_delta*y5+d*y5/I_r)-1/2*(y3-y2)*d*(y5/I_r)^2-1/2*(y3-y2)*d*(y6-y5/I_r)^2+1/2*d*((y3+y2)^2*(y5-d*y6)-2*(y3+1)*(y3-y2+1)*E_delta*y5*(y5-d*y6)+(y3-y2)^2*E_delta^2*y5^2))*h_r2+alpha*(y2*y5+(y3-y2)*(1-E_delta)*y5)*c_r+1/d*(y2^2*y5^2+2*y2*(y3-y2)*(1-E_delta)*y5+(y3-y2)^2*(1-E_delta)^2*y5^2)*h_c+((1-beta)*E_L*1/((1-l_3)*x1)*((u_1+E_L)^(1-l_3)-u_1^(1-l_3))+beta*E_L*y1*y6-(sigma/(1+l_1)*(2+l_1))*((1+trk1+y1*y6)^(2+l_1)-(1-trk1)^(2+l_1))-(sigma*y1*y6)/(1+l_1)*(1+trk1)^(1+l_1))*h_c)));

fp1(p,1)=f1(y1,y2,y3,y4,y5,y6);


f2=@(y1,y2,y3,y4,y5,y6) ((d/(y2*y5+(y3-y2)*(1-E_delta)*y5))*(a_re*(sigma^2/(1+2*l_1))*((1+trk1+y1*y6)^(1+2*l_1)-(1+trk1)^(1+2*l_1))-b_re*(sigma*(1+l_1))*((1+trk1+y1*y6)^(1+l_1)-(1+trk1)^(1+l_1))+c_re*y1*y6+a_me*(lambda^2/(1+2*l_2))*((1+tpk1+(y4-y1)*y6)^(1+2*l_2)-(1+tpk1)^(1+2*l_2))-b_me*(lambda*(1+l_2))*((1+tpk1+(y4-y1)*y6)^(1+l_2)-(1+tpk1)^(1+l_2))+c_me*(y4-y1)*y6+(2*y3*L_mr+2*L_rc+2*L_cm)*E*v_e+(1/2*y3*(y3+1)*E_delta*y5*y6+1/d*(y2^2*y5^2+2*y2*(y3-y2)*(1-E_delta)*y5+(y3-y2)^2*(1-E_delta)^2*y5^2))*w_e+(1-beta)*E_L*d_e));

fp2(p,1)=f2(y1,y2,y3,y4,y5,y6);

end

%%%Carbon cap and trade policy

for p=1:pop

if (fp2(p,1)>C_max)

 fp1(p,1)=fp1(p,1)+c_t*(fp2(p,1)-C_max);

else 
    fp1(p,1)=fp1(p,1)-c_t*(C_max-fp2(p,1));

end

end



Iter=zeros(maxgen,1);
%CVV=zeros(maxgen,1);

for gen=1:maxgen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% first pareto front
Sp=cell(pop,1);
np=zeros(pop,1);

F=cell(pop,1);

R=zeros(pop,1);

F{1,1}=[];

for p=1:pop
    Sp{p,1}=[]; np(p,1)=0; 

    for  q=1:pop 

       if fp1(p,1)<fp1(q,1) && fp2(p,1)<fp2(q,1)
            Sp{p,1}=[Sp{p,1},q];
      
        elseif fp1(q,1)<fp1(p,1) && fp2(q,1)<fp2(p,1)
         np(p,1)=np(p,1)+1;
      end

    end

if np(p,1)==0
    R(p,1)=1;

    F{1,1}=[F{1,1},p];
   
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Finding the other fronts

i=1;
l1=length(F{1,1});

while l1~=0

    Q=[];

for k=1:l1

    l2=length(Sp{F{i,1}(1,k),1});

    %if l2~=0
    
    for s=1:l2

        np(Sp{F{i,1}(1,k),1}(1,s),1)= np(Sp{F{i,1}(1,k),1}(1,s),1)-1; 

if np(Sp{F{i,1}(1,k),1}(1,s),1)==0

    R(Sp{F{i,1}(1,k),1}(1,s),1)=i+1;

    Q=[Q,Sp{F{i,1}(1,k),1}(1,s)];
end

    end

    %end

end

i=i+1;

F{i,1}=Q;

l1=length(F{i,1});

end



%scatter(fp1,fp2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finding  Diversity of each solution

div=zeros(pop,1);


%length of each fronts
lf=zeros(pop,1);
for i=1:pop
lf(i,1)=length(F{i,1});
end



i=1;

while lf(i,1)~=0

S1=[];
S2=[];
%indentifying the solution in each fronts
T1=[];

for k=1:lf(i,1)
    T1=[T1;F{i,1}(1,k)];
    z1=fp1(F{i,1}(1,k),1);
    S1=[S1;z1]; 
    z2=fp2(F{i,1}(1,k),1);
    S2=[S2;z2];
end

T2=T1;

%%sorting with respect to first objective
for k=1:lf(i,1)-1
    for k1=k+1:lf(i,1)
        if S1(k,1)>S1(k1,1)
            tmp=S1(k,1);
            S1(k,1)=S1(k1,1);
            S1(k1,1)=tmp;
            tmp=T1(k,1);
            T1(k,1)=T1(k1,1);
            T1(k1,1)=tmp;
      
        end
    end
end

%%sorting with respect to second objective

for k=1:lf(i,1)-1
    for k1=k+1:lf(i,1)
        if S2(k,1)>S2(k1,1)
            tmp=S2(k,1);
            S2(k,1)=S2(k1,1);
            S2(k1,1)=tmp;
            tmp=T2(k,1);
            T2(k,1)=T2(k1,1);
            T2(k1,1)=tmp;
        end
    end
end



div(T1(1,1),1)=inf; div(T1(lf(i,1),1),1)=inf;

for k=2:lf(i,1)-1
   div(T1(k,1),1)= div(T1(k,1),1)+((S1(k+1,1)-S1(k-1,1))/(S1(lf(i,1),1)-S1(1,1)));
   div(T2(k,1),1)= div(T2(k,1),1)+((S2(k+1,1)-S2(k-1,1))/(S2(lf(i,1),1)-S2(1,1)));
end


i=i+1;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%GA_crowded_binary_Tournament_selection
%id=GA_tournament_selection(pop,fp);

id=zeros(pop,1);

for p=1:pop

    q=p;

    while q==p
        r=rand;
        q=floor(r*pop);

        if q==0
            q=1;
           break
        end
    end


 if R(p,1)<R(q,1)

        id(p,1)=p;
       
 

 elseif R(q,1)<R(p,1)
      id(p,1)=q;
     
 

 elseif R(p,1)==R(q,1) && div(p,1)>div(q,1)

        id(p,1)=p;
  

elseif R(p,1)==R(q,1) && div(p,1)==div(q,1)

     %%%doubt randomly selected
     id(p,1)=p;

 else 
     id(p,1)=p;

end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GA_two_point_xover
% bc=GA_two_point_xover(pop,id,pc,tlen,bp);

%%selected two random (not same)from the mating pool for crossing
bc=zeros(pop,tlen);
k=0;
while k<pop
q=1;
p=1;
while (q==p)
    r=rand;
    p=floor(r*pop);
    if p==0
        p=1;
    end
    r=rand;
    q=floor(r*pop);
    if q==0
        q=1;
    end
end
p=id(p,1);
q=id(q,1);

r=rand;
if r<=pc

    %Generate two random crossing site
    cs1=1;
    cs2=1;
    while cs1==cs2
        r=rand;
        cs1=floor(r*(tlen-1));
        if cs1==0
            cs1=1;
        end
        r=rand;
        cs2=floor(r*(tlen-1));

        if cs2==0
            cs2=1;
        end

        if cs1>cs2
            tmp=cs1;
            cs1=cs2;
            cs2=tmp;
        end

    end

    %Generate two children by crossing the two parents
    for s=1:cs1
        bc(k+1,s)=bp(p,s);
        bc(k+2,s)=bp(q,s);
    end

    for s=cs1+1:cs2
        bc(k+1,s)=bp(q,s);
        bc(k+2,s)=bp(p,s);
    end

    for s=cs2+1:tlen
        bc(k+1,s)=bp(p,s);
        bc(k+2,s)=bp(q,s);
    end

else
    %%copy the two parent as two children
    for s=1:tlen
        bc(k+1,s)=bp(p,s);
        bc(k+2,s)=bp(q,s);
    end
end
k=k+2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GA_binary_mutation
% bc=GA_binary_mutation(pop,pm,tlen,bc);

xc=zeros(pop,nv);

for p=1:pop

    while (xc(p,1)>=xc(p,2)) || (xc(p,2)>=xc(p,4)) || (xc(p,4)>=xc(p,3))
        
        for s=1:tlen
        r=rand;
        if r<=pm
            if bc(p,s)==0
                bc(p,s)=1;
            else
                bc(p,s)=0;
            end
        end
        end
          c=0;
    for j=1:nv
        sum=0;
        for k=1:slen(j,1)
            sum=sum+bc(p,c+k)*(2^(slen(j,1)-k));
        end
        xc(p,j)=lmt(j,1)+sum*((lmt(j,2)-lmt(j,1))/(2^slen(j,1)-1));
        c=c+slen(j,1);
    end

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function_evalution 1
%function_evalution 2
% fc=user_funct(pop,xp);%function value for parents

fc1=zeros(pop,1);
fc2=zeros(pop,1);

for p=1:pop

y1=xc(p,1);
y2=xc(p,2);
y3=xc(p,3);
y4=xc(p,4);
y5=xc(p,5);
y6=xc(p,6);


A1=0;
S1=zeros;
for j=1:ceil(y1)
    S1(j)=((1+trk1+j*y6)^(2+l_1)-(1+trk1+(j-1)*y6)^(2+l_1));
A1=A1+S1(j);
end

A2=0;
S2=zeros;
for j=1:ceil(y1)
    S2(j)=((1+trk1+(j-1)*y6)^(1+l_1));
A2=A2+S2(j);
end



A3=0;
S3=zeros;
for j=1:ceil(y1)-1
    S3(j)=((1+trk1+j*y6)^(1+l_1)-(1+trk1)^(1+l_1));
A3=A3+S3(j);
end

A4=0;
S4=zeros;
for j=ceil(y2)+1:ceil(y4)-1
    S4(j)=((1+tpk1+(j-y1+1)*y6)^(2+l_2)-(1+tpk1+(j-y1)*y6)^(2+l_2));
A4=A4+S4(j);
end


f1=@(y1,y2,y3,y4,y5,y6)(((d/(y2*y5+(y3-y2)*(1-E_delta)*y5))*(s_m+(y3-y2)*y5*d_m+(y3-y2)*y6*gamma_m+y2*y6*gamma_r+(y3-y2)*E_delta*y5*c_m+y3*y5*s_mr+(sigma/((1+l_1)*(2+l_2))*A1-(sigma*y6)/(1+l_1)*A2+(sigma*y6)/(1+l_1)*A3+(y2-y1)*(sigma*y6)/(1+l_1)*((1+trk1+y1*y6)^(1+l_1)-(1+trk1)^(1+l_1))-1/2*y2*(y2-1)*y5*y6)*hm1+(lambda/((1+l_2)*(2+l_2))*((1+tpk1+(y2-y1+1)*y6)^(2+l_2)-(1+tpk1)^(2+l_2))-1/(1+l_2)*(y2-y1+1)*lambda*y6*(1+tpk1)^(1+l_2)+lambda/((1+l_2)*(2+l_2))*A4-1/2*(y4^2-y4-y2^2+y2)*(1+tpk1)^(1+l_2)+y6*(y3-y4)*lambda/(1+l_2)*((1+tpk1+(y4-y1)*y6)^(1+l_2)-(1+tpk1)^(1+l_2))-1/2*(y3-y2-1)*y5*y6)*hm2+A_r+y3*y5*w_r+(y3-y2)*y5*s_r_dash*1/I_r+1/2*y3*(y3+1)*E_delta*y5*y6*h_r1+(1/2*y2^2*(y5-d*y6)*y6+1/2*y2*y5*y6+y2*(y3-y2)*y6*(y5-d*y6)+1/2*(y3-y2+1)*(y3-y2)*y6*((1-E_delta)*y5-d*y6)+(y3-y2)*y5*y6-(y3-y2)*(y6-y5/I_r)*(E_delta*y5+d*y5/I_r)-1/2*(y3-y2)*d*(y5/I_r)^2-1/2*(y3-y2)*d*(y6-y5/I_r)^2+1/2*d*((y3+y2)^2*(y5-d*y6)-2*(y3+1)*(y3-y2+1)*E_delta*y5*(y5-d*y6)+(y3-y2)^2*E_delta^2*y5^2))*h_r2+alpha*(y2*y5+(y3-y2)*(1-E_delta)*y5)*c_r+1/d*(y2^2*y5^2+2*y2*(y3-y2)*(1-E_delta)*y5+(y3-y2)^2*(1-E_delta)^2*y5^2)*h_c+((1-beta)*E_L*1/((1-l_3)*x1)*((u_1+E_L)^(1-l_3)-u_1^(1-l_3))+beta*E_L*y1*y6-(sigma/(1+l_1)*(2+l_1))*((1+trk1+y1*y6)^(2+l_1)-(1-trk1)^(2+l_1))-(sigma*y1*y6)/(1+l_1)*(1+trk1)^(1+l_1))*h_c)));

fc1(p,1)=f1(y1,y2,y3,y4,y5,y6);


f2=@(y1,y2,y3,y4,y5,y6)((d/(y2*y5+(y3-y2)*(1-E_delta)*y5))*(a_re*(sigma^2/(1+2*l_1))*((1+trk1+y1*y6)^(1+2*l_1)-(1+trk1)^(1+2*l_1))-b_re*(sigma*(1+l_1))*((1+trk1+y1*y6)^(1+l_1)-(1+trk1)^(1+l_1))+c_re*y1*y6+a_me*(lambda^2/(1+2*l_2))*((1+tpk1+(y4-y1)*y6)^(1+2*l_2)-(1+tpk1)^(1+2*l_2))-b_me*(lambda*(1+l_2))*((1+tpk1+(y4-y1)*y6)^(1+l_2)-(1+tpk1)^(1+l_2))+c_me*(y4-y1)*y6+(2*y3*L_mr+2*L_rc+2*L_cm)*E*v_e+(1/2*y3*(y3+1)*E_delta*y5*y6+1/d*(y2^2*y5^2+2*y2*(y3-y2)*(1-E_delta)*y5+(y3-y2)^2*(1-E_delta)^2*y5^2))*w_e+(1-beta)*E_L*d_e));

fc2(p,1)=f2(y1,y2,y3,y4,y5,y6);

end



%%%Carbon cap and trade policy

for p=1:pop

if (fp2(p,1)>C_max)

 fp1(p,1)=fp1(p,1)+c_t*(fp2(p,1)-C_max);

else 
    fp1(p,1)=fp1(p,1)-c_t*(C_max-fp2(p,1));

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Elite_preservation
% [bp,xp,fp]=elite_preservation(pop,nv,tlen,bp,xp,fp,bc,xc,fc);


%%Form a combined population by mixing the parent and children populations


b=zeros(2*pop,tlen);
x=zeros(2*pop,nv);
f1=zeros(2*pop,1);
f2=zeros(2*pop,1);

for i=1:pop
    for s=1:tlen
        b(i,s)=bp(i,s);
    end
    for j=1:nv
        x(i,j)=xp(i,j);
    end

    f1(i,1)=fp1(i,1);
    f2(i,1)=fp2(i,1);

    for s=1:tlen
        b(pop+i,s)=bc(i,s);
    end
    for j=1:nv
        x(pop+i,j)=xc(i,j);
    end
    f1(pop+i,1)=fc1(i,1);
    f2(pop+i,1)=fc2(i,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v1=zeros(2*pop,1);
v2=zeros(2*pop,1);
v3=zeros(2*pop,1);
v4=zeros(2*pop,1);
v5=zeros(2*pop,1);
v6=zeros(2*pop,1);

for p=1:2*pop

y1=x(p,1);
y2=x(p,2);
y3=x(p,3);
y4=x(p,4);
y5=x(p,5);
y6=x(p,6);


%constraint 1 (equality)
g1=@(y1,y2,y3,y4,y5,y6) (sigma/(1+l_1))*((1+trk1+y1*y6)^(1+l_1)-(1+trk1)^(1+l_1))-y2*y5;


%constraint 2 (equality )
g2=@(y1,y2,y3,y4,y5,y6) (lambda/(1+l_2))*((1+tpk1+(y4-y1)*y6)^(1+l_1)-(1+tpk1)^(1+l_2))-(y3-y2)*y5;


%constraint 3 (less than type)
g3=@(y1,y2,y3,y4,y5,y6)  y5-(sigma/(1+l_1))*((1+trk1+y6)^(1+l_1)-(1+trk1)^(1+l_1));


%constraint 4 (less than type)
g4=@(y1,y2,y3,y4,y5,y6)  y6-(lambda/(1+l_2))*((1+tpk1+y6)^(1+l_1)-(1+tpk1)^(1+l_2));

%constraint 5 (less than  type)
g5=@(y1,y2,y3,y4,y5,y6)  d*y6-y5-E_delta*y5;


%constraint 6  (equality)
g6=@(y1,y2,y3,y4,y5,y6) (((1+l_1)*beta/sigma)*E_L+(1+trk1)^(1+l_1))^(1+l_1)-(1+trk1)-y1*y6;

v1(p,1)=g1(y1,y2,y3,y4,y5,y6);
v2(p,1)=g2(y1,y2,y3,y4,y5,y6);
v3(p,1)=g3(y1,y2,y3,y4,y5,y6);
v4(p,1)=g4(y1,y2,y3,y4,y5,y6);
v5(p,1)=g5(y1,y2,y3,y4,y5,y6);
v6(p,1)=g6(y1,y2,y3,y4,y5,y6);
end

fs1=zeros(2*pop,1);
fs2=zeros(2*pop,1);
% finding worst feasible solution f_max 
for p=1:2*pop
    if ((v1(p,1)==0) && (v2(p,1)==0) && (v3(p,1)<=0) && (v4(p,1)<=0)  && (v5(p,1)<=0) && (v6(p,1)==0))

        fs1(p,1)=f1(p,1);
        fs2(p,1)=f2(p,1);

    else
       fs1(p,1)=0; 
       fs2(p,1)=0; 

    end
end

%ascending
for p=1:2*pop-1
    for q=p+1:2*pop
        if fs1(p,1)>fs1(q,1)
            tmp=fs1(p,1);
            fs1(p,1)=fs1(q,1);
            fs1(q,1)=tmp;
        end
    end
end

%ascending
for p=1:2*pop-1
    for q=p+1:2*pop
        if fs2(p,1)>fs2(q,1)
            tmp=fs2(p,1);
            fs2(p,1)=fs2(q,1);
            fs2(q,1)=tmp;
        end
    end
end

%Total constraint violations
CV=zeros(2*pop,1);
for p=1:2*pop
     if (v1(p,1)~=0)
       CV(p,1)=abs(v1(p,1));
     else
       CV(p,1)=CV(p,1);
     end

     if (v2(p,1)~=0)
       CV(p,1)=CV(p,1)+abs(v2(p,1));
     else
        CV(p,1)=CV(p,1);
     end
     if (v3(p,1)>0)
     CV(p,1)=CV(p,1)+abs(v3(p,1));
     else
        CV(p,1)=CV(p,1);
     end

     if (v4(p,1)>0)
        CV(p,1)=CV(p,1)+abs(v4(p,1));
     else
       CV(p,1)=CV(p,1);
     end

     if (v5(p,1)>0)
        CV(p,1)=CV(p,1)+abs(v5(p,1));
     else
       CV(p,1)=CV(p,1);
     end

     if (v6(p,1)~=0)
        CV(p,1)=CV(p,1)+abs(v6(p,1));
     else
       CV(p,1)=CV(p,1);
     end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%check
%fitness value by Deb's approach

for p=1:2*pop
     if ((v1(p,1)==0) && (v2(p,1)==0) && (v3(p,1)<=0) && (v4(p,1)<=0)  && (v5(p,1)<=0) && (v6(p,1)==0))
       f1(p,1)=f1(p,1); 
       f2(p,1)=f2(p,1);

     else
         f1(p,1)=fs1(2*pop,1)+CV(p,1);
         f2(p,1)=fs2(2*pop,1)+CV(p,1);

     end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% first pareto front for combined population

Spc=cell(2*pop,1);
npc=zeros(2*pop,1);
Fc=cell(2*pop,1);
Rc=zeros(2*pop,1);
Fc{1,1}=[];

for p=1:2*pop

    Spc{p,1}=[]; npc(p,1)=0; 

    for  q=1:2*pop 
       if f1(p,1)<f1(q,1) && f2(p,1)<f2(q,1)
            Spc{p,1}=[Spc{p,1},q];
      
       elseif f1(q,1)<f1(p,1) && f2(q,1)<f2(p,1)
         npc(p,1)=npc(p,1)+1;
      end

    end

if npc(p,1)==0
    Rc(p,1)=1;

    Fc{1,1}=[Fc{1,1},p];
   
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Finding the other fronts

i=1;

l1c=length(Fc{1,1});

while l1c~=0

    Qc=[];

for k=1:l1c

    l2c=length(Spc{Fc{i,1}(1,k),1});

    %if l2~=0
    
    for s=1:l2c

        npc(Spc{Fc{i,1}(1,k),1}(1,s),1)= npc(Spc{Fc{i,1}(1,k),1}(1,s),1)-1; 

if npc(Spc{Fc{i,1}(1,k),1}(1,s),1)==0

    Rc(Spc{Fc{i,1}(1,k),1}(1,s),1)=i+1;

    Qc=[Qc,Spc{Fc{i,1}(1,k),1}(1,s)];
end

    end

    %end

end

i=i+1;

Fc{i,1}=Qc;

l1c=length(Fc{i,1});

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finding  Diversity of each solution

divc=zeros(2*pop,1);


%length of each fronts
lfc=zeros(2*pop,1);

for i=1:2*pop
lfc(i,1)=length(Fc{i,1});
end



i=1;

while lfc(i,1)~=0

S1c=[];
S2c=[];
%indentifying the solution in each fronts
T1c=[];

for k=1:lfc(i,1)
    T1c=[T1c;Fc{i,1}(1,k)];
    z1c=f1(Fc{i,1}(1,k),1);
    S1c=[S1c;z1c]; 
    z2c=f2(Fc{i,1}(1,k),1);
    S2c=[S2c;z2c];
end

T2c=T1c;

%%sorting with respect to first objective
for k=1:lfc(i,1)-1
    for k1=k+1:lfc(i,1)
        if S1c(k,1)>S1c(k1,1)
            tmp=S1c(k,1);
            S1c(k,1)=S1c(k1,1);
            S1c(k1,1)=tmp;
            tmp=T1c(k,1);
            T1c(k,1)=T1c(k1,1);
            T1c(k1,1)=tmp;
      
        end
    end
end

%%sorting with respect to second objective

for k=1:lfc(i,1)-1
    for k1=k+1:lfc(i,1)
        if S2c(k,1)>S2c(k1,1)
            tmp=S2c(k,1);
            S2c(k,1)=S2c(k1,1);
            S2c(k1,1)=tmp;
            tmp=T2c(k,1);
            T2c(k,1)=T2c(k1,1);
            T2c(k1,1)=tmp;
        end
    end
end



divc(T1c(1,1),1)=inf;

divc(T1c(lfc(i,1),1),1)=inf;

for k=2:lfc(i,1)-1
   divc(T1c(k,1),1)= divc(T1c(k,1),1)+((S1c(k+1,1)-S1c(k-1,1))/(S1c(lfc(i,1),1)-S1c(1,1)));
   divc(T2c(k,1),1)= divc(T2c(k,1),1)+((S2c(k+1,1)-S2c(k-1,1))/(S2c(lfc(i,1),1)-S2c(1,1)));
end


i=i+1;

end


U=zeros(pop,1);

p=0;
i=1;

while p<pop


lu=length(Fc{i,1});

for j=1:lu
    U(p+j,1)=Fc{i,1}(1,j);

end

i=i+1;
p=p+lu;

end 

V=zeros(pop,1);
for p=1:pop
    V(p,1)=U(p,1);
end

%Form the next population with the first 50% individuals of the combined population

for i=1:pop
    for s=1:tlen
        bp(i,s)=b(V(i,1),s);
    end
    for j=1:nv
        xp(i,j)=x(V(i,1),j);
    end
    fp1(i,1)=f1(V(i,1),1);
    fp2(i,1)=f2(V(i,1),1);

end

Iter(gen,1)=gen;
%CVV(gen,1)=CV(1,1);
fprintf('Iteration= %d\n',gen)

end 





























































 























































