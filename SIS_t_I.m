clear;
clc;

load('Degree_Distribution.mat');
M=max(uniqueDegree_1);%һ��facet������
N=max(uniqueDegree_2);%����facet������
M1=length(uniqueDegree_1);%һά��ĳ���
N1=length(uniqueDegree_2);%��ά��ĳ���
X0=zeros(1,2*M1*N1);
S0=1900; I0=100;%��ʼʱ�������׸��ߺ����и�Ⱦ�ߵ��ܶ�
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i����һά��
    for j=1:N1 %j�����ά��
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end
[T11,X1]=ode45(@SIS, [0 300] ,X0);
%��ÿ��ʱ�̵�����I������
I1=[];
for i=1:length(T11)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum+X1(i,j);
    end
    I1(i)=sum;
end
save I0=100;
%%
clear;
clc;

load('Degree_Distribution.mat');
M=max(uniqueDegree_1);%һ��facet������
N=max(uniqueDegree_2);%����facet������
M1=length(uniqueDegree_1);%һά��ĳ���
N1=length(uniqueDegree_2);%��ά��ĳ���
X0=zeros(1,2*M1*N1);
S0=1800; I0=200;%��ʼʱ�������׸��ߺ����и�Ⱦ�ߵ��ܶ�
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i����һά��
    for j=1:N1 %j�����ά��
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end
[T12,X1]=ode45(@SIS, [0 300] ,X0);
%��ÿ��ʱ�̵�����I������
I2=[];
for i=1:length(T12)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum+X1(i,j);
    end
    I2(i)=sum;
end
save I0=200;
%%
clear;
clc;

load('Degree_Distribution.mat');
M=max(uniqueDegree_1);%һ��facet������
N=max(uniqueDegree_2);%����facet������
M1=length(uniqueDegree_1);%һά��ĳ���
N1=length(uniqueDegree_2);%��ά��ĳ���
X0=zeros(1,2*M1*N1);
S0=1600; I0=400;%��ʼʱ�������׸��ߺ����и�Ⱦ�ߵ��ܶ�
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i����һά��
    for j=1:N1 %j�����ά��
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end
[T13,X1]=ode45(@SIS, [0 300] ,X0);
%��ÿ��ʱ�̵�����I������
I3=[];
for i=1:length(T13)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum+X1(i,j);
    end
    I3(i)=sum;
end
save I0=400;
%%
clear;
clc;

load('Degree_Distribution.mat');
M=max(uniqueDegree_1);%һ��facet������
N=max(uniqueDegree_2);%����facet������
M1=length(uniqueDegree_1);%һά��ĳ���
N1=length(uniqueDegree_2);%��ά��ĳ���
X0=zeros(1,2*M1*N1);
S0=1400; I0=600;%��ʼʱ�������׸��ߺ����и�Ⱦ�ߵ��ܶ�
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i����һά��
    for j=1:N1 %j�����ά��
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end
[T14,X1]=ode45(@SIS, [0 300] ,X0);
%��ÿ��ʱ�̵�����I������
I4=[];
for i=1:length(T14)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum+X1(i,j);
    end
    I4(i)=sum;
end
save I0=600;
%%
clear;
clc;

load('Degree_Distribution.mat');
M=max(uniqueDegree_1);%һ��facet������
N=max(uniqueDegree_2);%����facet������
M1=length(uniqueDegree_1);%һά��ĳ���
N1=length(uniqueDegree_2);%��ά��ĳ���
X0=zeros(1,2*M1*N1);
S0=1200; I0=800;%��ʼʱ�������׸��ߺ����и�Ⱦ�ߵ��ܶ�
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i����һά��
    for j=1:N1 %j�����ά��
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end
[T15,X1]=ode45(@SIS, [0 300] ,X0);
%��ÿ��ʱ�̵�����I������
I5=[];
for i=1:length(T15)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum+X1(i,j);
    end
    I5(i)=sum;
end
save I0=800;
%%
clear;
clc;

load('Degree_Distribution.mat');
M=max(uniqueDegree_1);%һ��facet������
N=max(uniqueDegree_2);%����facet������
M1=length(uniqueDegree_1);%һά��ĳ���
N1=length(uniqueDegree_2);%��ά��ĳ���
X0=zeros(1,2*M1*N1);
S0=1000; I0=1000;%��ʼʱ�������׸��ߺ����и�Ⱦ�ߵ��ܶ�
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i����һά��
    for j=1:N1 %j�����ά��
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end
[T16,X1]=ode45(@SIS, [0 300] ,X0);
%��ÿ��ʱ�̵�����I������
I6=[];
for i=1:length(T16)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum+X1(i,j);
    end
    I6(i)=sum;
end
save I0=1000;
%%
clear;
clc;

load('Degree_Distribution.mat');
M=max(uniqueDegree_1);%һ��facet������
N=max(uniqueDegree_2);%����facet������
M1=length(uniqueDegree_1);%һά��ĳ���
N1=length(uniqueDegree_2);%��ά��ĳ���
X0=zeros(1,2*M1*N1);
S0=800; I0=1200;%��ʼʱ�������׸��ߺ����и�Ⱦ�ߵ��ܶ�
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i����һά��
    for j=1:N1 %j�����ά��
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end
[T17,X1]=ode45(@SIS, [0 300] ,X0);
%��ÿ��ʱ�̵�����I������
I7=[];
for i=1:length(T17)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum+X1(i,j);
    end
    I7(i)=sum;
end
save I0=1200;
%%
clear;
clc;

load('Degree_Distribution.mat');
M=max(uniqueDegree_1);%һ��facet������
N=max(uniqueDegree_2);%����facet������
M1=length(uniqueDegree_1);%һά��ĳ���
N1=length(uniqueDegree_2);%��ά��ĳ���
X0=zeros(1,2*M1*N1);
S0=600; I0=1400;%��ʼʱ�������׸��ߺ����и�Ⱦ�ߵ��ܶ�
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i����һά��
    for j=1:N1 %j�����ά��
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end
[T18,X1]=ode45(@SIS, [0 300] ,X0);
%��ÿ��ʱ�̵�����I������
I8=[];
for i=1:length(T18)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum+X1(i,j);
    end
    I8(i)=sum;
end
save I0=1400;
%%
clear;
clc;

load('Degree_Distribution.mat');
M=max(uniqueDegree_1);%һ��facet������
N=max(uniqueDegree_2);%����facet������
M1=length(uniqueDegree_1);%һά��ĳ���
N1=length(uniqueDegree_2);%��ά��ĳ���
X0=zeros(1,2*M1*N1);
S0=400; I0=1600;%��ʼʱ�������׸��ߺ����и�Ⱦ�ߵ��ܶ�
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i����һά��
    for j=1:N1 %j�����ά��
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end
[T19,X1]=ode45(@SIS, [0 300] ,X0);
%��ÿ��ʱ�̵�����I������
I9=[];
for i=1:length(T19)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum+X1(i,j);
    end
    I9(i)=sum;
end
save I0=1600;
%%
clear;
clc;

load('Degree_Distribution.mat');
M=max(uniqueDegree_1);%һ��facet������
N=max(uniqueDegree_2);%����facet������
M1=length(uniqueDegree_1);%һά��ĳ���
N1=length(uniqueDegree_2);%��ά��ĳ���
X0=zeros(1,2*M1*N1);
S0=200; I0=1800;%��ʼʱ�������׸��ߺ����и�Ⱦ�ߵ��ܶ�
%X(0)=[S00,I00,S01,I01,S02,I02,...,S0N1,I0N1,...,SM10,IM10,SM11,IM11,...,SM1N1,IM1N1]
for i=1:M1 %i����һά��
    for j=1:N1 %j�����ά��
        X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
        X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
    end
end
[T110,X1]=ode45(@SIS, [0 300] ,X0);
%��ÿ��ʱ�̵�����I������
I10=[];
for i=1:length(T110)
    sum=0;
    for j=2:2:length(X1(1,:))
        sum = sum+X1(i,j);
    end
    I10(i)=sum;
end
save I0=1800;

%%
load('I0=100.mat');
load('I0=200.mat');
load('I0=400.mat');
load('I0=600.mat');
load('I0=800.mat');
load('I0=1000.mat');
load('I0=1200.mat');
load('I0=1400.mat');
load('I0=1600.mat');
load('I0=1800.mat');
plot(T11,I1,'linewidth',2);
hold on;
plot(T12,I2,'linewidth',2);
hold on;
plot(T13,I3,'linewidth',2);
hold on;
plot(T14,I4,'linewidth',2);
hold on;
plot(T15,I5,'linewidth',2);
hold on;
plot(T16,I6,'linewidth',2);
hold on;
plot(T17,I7,'linewidth',2);
hold on;
plot(T18,I8,'linewidth',2);
hold on;
plot(T19,I9,'linewidth',2);
hold on;
plot(T110,I10,'linewidth',2);
hold on;

xlabel('t');
ylabel('I');
set(gca,'color','none');
hold on;

