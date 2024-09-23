%ctrl+R一次注释多行
%ctrl+T取消注释
%ctrl+C终止
%求解微分方程组的一种方法，4阶五级Runge-Kutta算法
%[t,x]=ode45(Fun,tspan,x0,options,pars)
%t是列向量，X为n×n矩阵，它的每一列就是其中一个状态变量随t的变化值
%fun为要求解的方程组，形式必须为一阶微分方程组；@代表方程组右边的部分，并且是列向量
%tspan为要求解的t的区间，x0为初始值，options是一些选项
clear;
clc;
global beta1 beta2 gamma %定义全局变量（必须，否则会出错）
gamma=0.15;beta2=0.0625;
load('Degree_Distribution.mat');

%%
Beta11=[0:0.001:0.03];%给定beta1要取的范围
Infect_11=zeros(1,length(Beta11));%每个beta1都会对应一个I，所以需要给I申请一个空间来存放
%通过循环将每个beta1对应的I算出来并放在Infect中
for i=1:length(Beta11)
    beta1=Beta11(i);
    M1=length(uniqueDegree_1);%一维面的长度
    N1=length(uniqueDegree_2);%二维面的长度
    X0=zeros(1,2*M1*N1);
    S0=1950; I0=50;%初始时刻所有易感者和所有感染者的数量
    for k=1:M1 %i代表一维面
        for j=1:N1 %j代表二维面
            X0(2*(k-1)*N1+2*j-1)=S0*Com_num(j,k);
            X0(2*(k-1)*N1+2*j)=I0*Com_num(j,k);
        end
    end
    [t,x]=ode45(@SIS_fun,[0 225],X0,[]);
    sum=0;
    for j=2:2:length(x(1,:))
        sum = sum+x(end,j);
    end
    Infect_11(i)=sum;
    fprintf('beta1=%f时完成！\n',Beta11(i));   
end
h1=plot(Beta11,Infect_11,'-b','linewidth',2);
hold on
%%
Beta1=[0:0.001:0.03];%给定beta1要取的范围
Infect=zeros(1,length(Beta1));%每个beta1都会对应一个I，所以需要给I申请一个空间来存放
%通过循环将每个beta1对应的I算出来并放在Infect中
for i=1:length(Beta1)
    beta1=Beta1(i);
    M1=length(uniqueDegree_1);%一维面的长度
    N1=length(uniqueDegree_2);%二维面的长度
    X0=zeros(1,2*M1*N1);
    S0=1500; I0=500;%初始时刻所有易感者和所有感染者的密度
    for k=1:M1 %i代表一维面
        for j=1:N1 %j代表二维面
            X0(2*(k-1)*N1+2*j-1)=S0*Com_num(j,k);
            X0(2*(k-1)*N1+2*j)=I0*Com_num(j,k);
        end
    end
    [t,x]=ode45(@SIS_fun,[0 225],X0,[]);
    sum=0;
    for j=2:2:length(x(1,:))
        sum = sum+x(end,j);
    end
    Infect(i)=sum;
    fprintf('beta1=%f时完成！\n',Beta1(i));
end
h2=plot(Beta1,Infect,'-.r','linewidth',1);
hold on

%%
%图例、"\phi"代表输出希腊字母，_代表输出下标，\infty代表输出无穷
legend([h1,h2],'I=50','I=500','Location','NorthWest');
xlabel('\beta_1');
ylabel('I');
%设置背景透明
set(gca,'color','none');
save data_beta1_N_2000_I_50_500_U





