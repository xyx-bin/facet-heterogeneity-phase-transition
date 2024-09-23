clear;
clc;

global beta1;
beta1 = 0.003;

%% 
load('facet1=Zipf_facet2=poissrnd_2.mat');
M=max(uniqueDegree_1_ult);%一阶facet的最大度
N=max(uniqueDegree_2_ult);%二阶facet的最大度
M1=length(uniqueDegree_1_ult);%一维面的长度
N1=length(uniqueDegree_2_ult);%二维面的长度
X0=zeros(1,2*M1*N1);
S0=1200; I0=800;%初始时刻所有易感者和所有感染者的密度
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i代表一维面
    for j=1:N1 %j代表二维面
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end

[T1,X1]=ode45(@distribution_fun2, [0:5:200] ,X0);
%把每个时刻的所有I加起来
I1=[];
for i=1:length(T1)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum + X1(i,j);
    end
    I1(i)=sum;
end
plot(T1,I1,'-ob','linewidth',1);
hold on;
h7=plot(T1,I1,'b','linewidth',1);
hold on;

%%
S0=1900; I0=100;%初始时刻所有易感者和所有感染者的密度
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i代表一维面
    for j=1:N1 %j代表二维面
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end

[T2,X2]=ode45(@distribution_fun2, [0:5:200] ,X0);
%把每个时刻的所有I加起来
I2=[];
for i=1:length(T2)
    sum=0;
    for j=2:2:length(X2(1,:))
        sum = sum + X2(i,j);
    end
    I2(i)=sum;
end
h1=plot(T2,I2,'-ob','linewidth',1);
hold on;

%%
load('facet1=Zipf_facet2=poissrnd_3.mat');
M=max(uniqueDegree_1_ult);%一阶facet的最大度
N=max(uniqueDegree_2_ult);%二阶facet的最大度
M1=length(uniqueDegree_1_ult);%一维面的长度
N1=length(uniqueDegree_2_ult);%二维面的长度
X0=zeros(1,2*M1*N1);
S0=1200; I0=800;
for i=1:M1
    for j=1:N1 
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end

[T3,X3]=ode45(@distribution_fun3, [0:5:200] ,X0);
I3=[];
for i=1:length(T3)
    sum=0;
    for j=2:2:length(X3(1,:))
        sum = sum + X3(i,j);
    end
    I3(i)=sum;
end
plot(T3,I3,'-*b','linewidth',1);
hold on;

%%
S0=1900; I0=100;
for i=1:M1 
    for j=1:N1 
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end

[T4,X4]=ode45(@distribution_fun3, [0:5:200] ,X0);
I4=[];
for i=1:length(T4)
    sum=0;
    for j=2:2:length(X4(1,:))
        sum = sum + X4(i,j);
    end
    I4(i)=sum;
end
h2=plot(T4,I4,'-*b','linewidth',1);
hold on;

%%
load('facet1=Zipf_facet2=poissrnd_4.mat');
M=max(uniqueDegree_1_ult);%一阶facet的最大度
N=max(uniqueDegree_2_ult);%二阶facet的最大度
M1=length(uniqueDegree_1_ult);%一维面的长度
N1=length(uniqueDegree_2_ult);%二维面的长度
X0=zeros(1,2*M1*N1);
S0=1200; I0=800;
for i=1:M1
    for j=1:N1 
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end

[T5,X5]=ode45(@distribution_fun4, [0:5:200] ,X0);
I5=[];
for i=1:length(T5)
    sum=0;
    for j=2:2:length(X5(1,:))
        sum = sum + X5(i,j);
    end
    I5(i)=sum;
end
plot(T5,I5,'-sb','linewidth',1);
hold on;

%%
S0=1900; I0=100;
for i=1:M1 
    for j=1:N1 
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end

[T6,X6]=ode45(@distribution_fun4, [0:5:200] ,X0);
I6=[];
for i=1:length(T6)
    sum=0;
    for j=2:2:length(X6(1,:))
        sum = sum + X6(i,j);
    end
    I6(i)=sum;
end
h3=plot(T6,I6,'-sb','linewidth',1);
hold on;

%% 
%Base Model的t_I
global beta2 gamma k1 k2;
beta2=0.0625;
gamma=0.15;
k1=20;
k2=6;

%%
X0=0.05;
[T7,X7]=ode45(@myfun, [0:5:200] ,X0);
%把每个时刻的所有I加起
for i=1:length(X7)
    X7(i)=X7(i)*2000;
end
h4=plot(T7,X7,'-or','linewidth',1);
hold on;
%%
X0=0.4;
[T8,X8]=ode45(@myfun, [0:5:200] ,X0);
%把每个时刻的所有I加起
for i=1:length(X8)
    X8(i)=X8(i)*2000;
end
plot(T8,X8,'-or','linewidth',1);
hold on;

%%
X0=0.05;
[T9,X9]=ode45(@myfun, [0:8:200] ,X0);
for i=1:length(X9)
    X9(i)=X9(i)*2000;
end
h5=plot(T9,X9,'-*r','linewidth',1);
hold on;
%%
X0=0.4;
[T10,X10]=ode45(@myfun,[0:8:200] ,X0);
for i=1:length(X10)
    X10(i)=X10(i)*2000;
end
plot(T10,X10,'-*r','linewidth',1);
hold on;

%%
X0=0.05;
[T11,X11]=ode45(@myfun, [0:2:200] ,X0);
for i=1:length(X11)
    X11(i)=X11(i)*2000;
end
h6=plot(T11,X11,'-sr','linewidth',1);
hold on;
h8=plot(T11,X11,'r','linewidth',1);
hold on;
%%
X0=0.4;
[T12,X12]=ode45(@myfun, [0:2:200] ,X0);
for i=1:length(X12)
    X12(i)=X12(i)*2000;
end
plot(T12,X12,'-sr','linewidth',1);
hold on;

legend([h7,h8,h1,h2,h3],'General Model','Base Model','<k_1^2>/<k_1>=22.87','<k_1^2>/<k_1>=27.59','<k_1^2>/<k_1>=42.88','Location','NorthWest');

xlabel('t');
ylabel('I');
set(gca,'color','none');
hold on;
