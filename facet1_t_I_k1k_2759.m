clear;
clc;

global beta1;
beta1 = 0.004;
load('facet1=Zipf_facet2=poissrnd_3.mat');
M=max(uniqueDegree_1_ult);%一阶facet的最大度
N=max(uniqueDegree_2_ult);%二阶facet的最大度
M1=length(uniqueDegree_1_ult);%一维面的长度
N1=length(uniqueDegree_2_ult);%二维面的长度
X0=zeros(1,2*M1*N1);

%%
S0=200; I0=1800;%初始时刻所有易感者和所有感染者的密度
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i代表一维面
    for j=1:N1 %j代表二维面
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end

[T1,X1]=ode45(@distribution_fun3, [0 200] ,X0);
%把每个时刻的所有I加起来
I1=[];
for i=1:length(T1)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum + X1(i,j);
    end
    I1(i)=sum;
end
plot(T1,I1,'linewidth',2);
hold on;
%%
S0=400; I0=1600;%初始时刻所有易感者和所有感染者的密度
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i代表一维面
    for j=1:N1 %j代表二维面
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end

[T2,X1]=ode45(@distribution_fun3, [0 200] ,X0);
%把每个时刻的所有I加起来
I2=[];
for i=1:length(T2)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum + X1(i,j);
    end
    I2(i)=sum;
end
plot(T2,I2,'linewidth',2);
hold on;
%%
S0=600; I0=1400;%初始时刻所有易感者和所有感染者的密度
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i代表一维面
    for j=1:N1 %j代表二维面
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end

[T3,X1]=ode45(@distribution_fun3, [0 200] ,X0);
%把每个时刻的所有I加起来
I3=[];
for i=1:length(T3)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum + X1(i,j);
    end
    I3(i)=sum;
end
plot(T3,I3,'linewidth',2);
hold on;
%%
S0=800; I0=1200;%初始时刻所有易感者和所有感染者的密度
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i代表一维面
    for j=1:N1 %j代表二维面
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end

[T4,X1]=ode45(@distribution_fun3, [0 200] ,X0);
%把每个时刻的所有I加起来
I4=[];
for i=1:length(T4)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum + X1(i,j);
    end
    I4(i)=sum;
end
plot(T4,I4,'linewidth',2);
hold on;
%%
S0=1000; I0=1000;%初始时刻所有易感者和所有感染者的密度
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i代表一维面
    for j=1:N1 %j代表二维面
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end

[T5,X1]=ode45(@distribution_fun3, [0 200] ,X0);
%把每个时刻的所有I加起来
I5=[];
for i=1:length(T5)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum + X1(i,j);
    end
    I5(i)=sum;
end
plot(T5,I5,'linewidth',2);
hold on;
%%
S0=1200; I0=800;%初始时刻所有易感者和所有感染者的密度
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i代表一维面
    for j=1:N1 %j代表二维面
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end

[T6,X1]=ode45(@distribution_fun3, [0 200] ,X0);
%把每个时刻的所有I加起来
I6=[];
for i=1:length(T6)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum + X1(i,j);
    end
    I6(i)=sum;
end
plot(T6,I6,'linewidth',2);
hold on;
%%
S0=1400; I0=600;%初始时刻所有易感者和所有感染者的密度
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i代表一维面
    for j=1:N1 %j代表二维面
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end

[T7,X1]=ode45(@distribution_fun3, [0 200] ,X0);
%把每个时刻的所有I加起来
I7=[];
for i=1:length(T7)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum + X1(i,j);
    end
    I7(i)=sum;
end
plot(T7,I7,'linewidth',2);
hold on;
%%
S0=1600; I0=400;%初始时刻所有易感者和所有感染者的密度
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i代表一维面
    for j=1:N1 %j代表二维面
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end

[T8,X1]=ode45(@distribution_fun3, [0 200] ,X0);
%把每个时刻的所有I加起来
I8=[];
for i=1:length(T8)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum + X1(i,j);
    end
    I8(i)=sum;
end
plot(T8,I8,'linewidth',2);
hold on;
%%
S0=1800; I0=200;%初始时刻所有易感者和所有感染者的密度
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i代表一维面
    for j=1:N1 %j代表二维面
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end

[T9,X1]=ode45(@distribution_fun3, [0 200] ,X0);
%把每个时刻的所有I加起来
I9=[];
for i=1:length(T9)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum + X1(i,j);
    end
    I9(i)=sum;
end
plot(T9,I9,'linewidth',2);
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

[T10,X1]=ode45(@distribution_fun3, [0 200] ,X0);
%把每个时刻的所有I加起来
I10=[];
for i=1:length(T10)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum + X1(i,j);
    end
    I10(i)=sum;
end
plot(T10,I10,'linewidth',2);
hold on;
%%
S0=1950; I0=50;%初始时刻所有易感者和所有感染者的密度
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i代表一维面
    for j=1:N1 %j代表二维面
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end

[T11,X1]=ode45(@distribution_fun3, [0 200] ,X0);
%把每个时刻的所有I加起来
I11=[];
for i=1:length(T11)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum + X1(i,j);
    end
    I11(i)=sum;
end
plot(T11,I11,'linewidth',2);
hold on;


xlabel('t');
ylabel('I');
set(gca,'color','none');
hold on;
