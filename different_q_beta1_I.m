clc;
clear;
global beta1 beta2 k1 k2 gamma q theta2 N average_2 %定义全局变量，C为3*三角形的数量=网络中的总边数
gamma=0.6;beta2=0.55;N=500;average_2=3; %average_2=3为二维面的平均度，k1、k2为节点连接三角形的个数
%D1=10;
%D2=10;
%beta1=0.01;
%k1=4;
%q=0.1;
%k2=round((average_2-k1*q)/(1-q)); % round四舍五入取整
%%
q=0.1;
k1=3;
k2=round((average_2-k1*q)/(1-q)); % round四舍五入取整
Infect1=[];
Beta1=[0:0.0025:0.3];
for i=1:length(Beta1)
    beta1=Beta1(i);
    time=0:0.1:100;
    x0=[]; %time初始值（即0的时候）对应的dx=[dE;dI]中E的初始值和I的初始值
    S0=480;I0=20;
    x0(1)=S0*q;
    x0(2)=I0*q;
    x0(3)=S0*(1-q);
    x0(4)=I0*(1-q);
    [t,x]=ode45(@simplify_fun1,time,x0,[]);
    Infect1(i)=x(end,2)+x(end,4);
    Beta1(i)
end

h1=plot(Beta1,Infect1,'-k','linewidth',1);
hold on
%%
q=0.1;
k1=4;
k2=round((average_2-k1*q)/(1-q)); % round四舍五入取整
Infect2=[];
Beta1=[0:0.0025:0.3];
for i=1:length(Beta1)
    beta1=Beta1(i);
    time=0:1:100;
    x0=[]; %time初始值（即0的时候）对应的dx=[dE;dI]中E的初始值和I的初始值
    S0=480;I0=20;
    x0(1)=S0*q;
    x0(2)=I0*q;
    x0(3)=S0*(1-q);
    x0(4)=I0*(1-q);
    [t,x]=ode45(@simplify_fun1,time,x0,[]);
    Infect2(i)=x(end,2)+x(end,4);
    Beta1(i)
end

h2=plot(Beta1,Infect2,'-r','linewidth',1);
hold on
xlabel('\beta_1');
ylabel('I');
%设置背景透明
set(gca,'color','none');

%%
q=0.1;
k1=5;
k2=round((average_2-k1*q)/(1-q)); % round四舍五入取整
Infect3=[];
Beta1=[0:0.0025:0.3];
for i=1:length(Beta1)
    beta1=Beta1(i);
    time=0:0.1:100;
    x0=[]; %time初始值（即0的时候）对应的dx=[dE;dI]中E的初始值和I的初始值
    S0=480;I0=20;
    x0(1)=S0*q;
    x0(2)=I0*q;
    x0(3)=S0*(1-q);
    x0(4)=I0*(1-q);
    [t,x]=ode45(@simplify_fun1,time,x0,[]);
    Infect3(i)=x(end,2)+x(end,4);
    Beta1(i)
end

h3=plot(Beta1,Infect3,'-g','linewidth',1);
hold on

%%
q=0.1;
k1=6;
k2=round((average_2-k1*q)/(1-q)); % round四舍五入取整
Infect4=[];
Beta1=[0:0.0025:0.3];
for i=1:length(Beta1)
    beta1=Beta1(i);
    time=0:0.1:100;
    x0=[]; %time初始值（即0的时候）对应的dx=[dE;dI]中E的初始值和I的初始值
    S0=480;I0=20;
    x0(1)=S0*q;
    x0(2)=I0*q;
    x0(3)=S0*(1-q);
    x0(4)=I0*(1-q);
    [t,x]=ode45(@simplify_fun1,time,x0,[]);
    Infect4(i)=x(end,2)+x(end,4);
    Beta1(i)
end

h4=plot(Beta1,Infect4,'-b','linewidth',1);
hold on


%%=============================================================================================
%%
q=0.1;
k1=3;
k2=round((average_2-k1*q)/(1-q)); % round四舍五入取整
Infect5=[];
Beta1=[0:0.0025:0.3];
for i=1:length(Beta1)
    beta1=Beta1(i);
    time=0:0.1:100;
    x0=[]; %time初始值（即0的时候）对应的dx=[dE;dI]中E的初始值和I的初始值
    S0=300;I0=200;
    x0(1)=S0*q;
    x0(2)=I0*q;
    x0(3)=S0*(1-q);
    x0(4)=I0*(1-q);
    [t,x]=ode45(@simplify_fun1,time,x0,[]);
    Infect5(i)=x(end,2)+x(end,4);
    Beta1(i)
end

h5=plot(Beta1,Infect5,'-k','linewidth',1);
hold on
%%
q=0.1;
k1=4;
k2=round((average_2-k1*q)/(1-q)); % round四舍五入取整
Infect6=[];
Beta1=[0:0.0025:0.3];
for i=1:length(Beta1)
    beta1=Beta1(i);
    time=0:1:100;
    x0=[]; %time初始值（即0的时候）对应的dx=[dE;dI]中E的初始值和I的初始值
    S0=300;I0=200;
    x0(1)=S0*q;
    x0(2)=I0*q;
    x0(3)=S0*(1-q);
    x0(4)=I0*(1-q);
    [t,x]=ode45(@simplify_fun1,time,x0,[]);
    Infect6(i)=x(end,2)+x(end,4);
    Beta1(i)
end

h6=plot(Beta1,Infect6,'-r','linewidth',1);
hold on

%%
q=0.1;
k1=5;
k2=round((average_2-k1*q)/(1-q)); % round四舍五入取整
Infect7=[];
Beta1=[0:0.0025:0.3];
for i=1:length(Beta1)
    beta1=Beta1(i);
    time=0:0.1:100;
    x0=[]; %time初始值（即0的时候）对应的dx=[dE;dI]中E的初始值和I的初始值
    S0=300;I0=200;
    x0(1)=S0*q;
    x0(2)=I0*q;
    x0(3)=S0*(1-q);
    x0(4)=I0*(1-q);
    [t,x]=ode45(@simplify_fun1,time,x0,[]);
    Infect7(i)=x(end,2)+x(end,4);
    Beta1(i)
end

h7=plot(Beta1,Infect7,'-g','linewidth',1);
hold on
%%
q=0.1;
k1=6;
k2=round((average_2-k1*q)/(1-q)); % round四舍五入取整
Infect8=[];
Beta1=[0:0.0025:0.3];
for i=1:length(Beta1)
    beta1=Beta1(i);
    time=0:0.1:100;
    x0=[]; %time初始值（即0的时候）对应的dx=[dE;dI]中E的初始值和I的初始值
    S0=300;I0=200;
    x0(1)=S0*q;
    x0(2)=I0*q;
    x0(3)=S0*(1-q);
    x0(4)=I0*(1-q);
    [t,x]=ode45(@simplify_fun1,time,x0,[]);
    Infect8(i)=x(end,2)+x(end,4);
    Beta1(i)
end

h8=plot(Beta1,Infect8,'-b','linewidth',1);
hold on

%%
legend([h1,h2,h7,h8],'k^2/k=3','k^2/k=3.13','k^2/k=3.31','k^2/k=3.55','Location','NorthWest');


%%
%function [y1,...,yN] = myfun(x1,...,xM) .名为 myfun 的函数，
%该函数接受输入 x1,...,xM 并返回输出 y1,...,yN
function dx=simplify_fun1(t,x)
%修改变量为全局变量
global beta1 beta2 k1 k2 gamma q N theta2 average_2
dx=zeros(4,1);
theta2=(k1*x(2)+k2*x(4))/(average_2*N);
dx(1)=-2*beta1*k1*x(1)*theta2-beta2*k1*x(1)*theta2^2+gamma*x(2);
dx(2)=2*beta1*k1*x(1)*theta2+beta2*k1*x(1)*theta2^2-gamma*x(2);
dx(3)=-2*beta1*k2*x(3)*theta2-beta2*k2*x(3)*theta2^2+gamma*x(4);
dx(4)=2*beta1*k2*x(3)*theta2+beta2*k2*x(3)*theta2^2-gamma*x(4);
end