clc;
clear;
global beta1 beta2 k1 k2 gamma q theta2 N average_2%定义全局变量
gamma=0.6;beta2=0.55;N=500;average_2=3;

% 调这些参数
beta1=0.08;
q=0.1;
k1=6;
k2=round((average_2-k1*q)/(1-q)); % round四舍五入取整
%%
Infect1=[];
time=0:0.5:150;
x0=[]; %time初始值（即0的时候）对应的dx=[dE;dI]中E的初始值和I的初始值
S0=480;I0=20;
x0(1)=S0*q;
x0(2)=I0*q;
x0(3)=S0*(1-q);
x0(4)=I0*(1-q);
[t,x]=ode45(@simplify_fun1,time,x0,[]);
for i=1:length(t)
    sum=0;
    for j=2:2:length(x(1,:))
        sum = sum+x(i,j);
    end
    Infect1(i)=sum;
end
h1=plot(t,Infect1,'linewidth',2);
hold on
%%
Infect2=[];
time=0:0.5:150;
x0=[]; %time初始值（即0的时候）对应的dx=[dE;dI]中E的初始值和I的初始值
S0=450;I0=50;
x0(1)=S0*q;
x0(2)=I0*q;
x0(3)=S0*(1-q);
x0(4)=I0*(1-q);
[t,x]=ode45(@simplify_fun1,time,x0,[]);
for i=1:length(t)
    sum=0;
    for j=2:2:length(x(1,:))
        sum = sum+x(i,j);
    end
    Infect2(i)=sum;
end
h2=plot(t,Infect2,'linewidth',2);
hold on
%%
Infect3=[];
time=0:0.5:150;
x0=[]; %time初始值（即0的时候）对应的dx=[dE;dI]中E的初始值和I的初始值
S0=400;I0=100;
x0(1)=S0*q;
x0(2)=I0*q;
x0(3)=S0*(1-q);
x0(4)=I0*(1-q);
[t,x]=ode45(@simplify_fun1,time,x0,[]);
for i=1:length(t)
    sum=0;
    for j=2:2:length(x(1,:))
        sum = sum+x(i,j);
    end
    Infect3(i)=sum;
end
h3=plot(t,Infect3,'linewidth',2);
hold on
%%
Infect4=[];
time=0:0.5:150;
x0=[]; %time初始值（即0的时候）对应的dx=[dE;dI]中E的初始值和I的初始值
S0=350;I0=150;
x0(1)=S0*q;
x0(2)=I0*q;
x0(3)=S0*(1-q);
x0(4)=I0*(1-q);
[t,x]=ode45(@simplify_fun1,time,x0,[]);
for i=1:length(t)
    sum=0;
    for j=2:2:length(x(1,:))
        sum = sum+x(i,j);
    end
    Infect4(i)=sum;
end
h4=plot(t,Infect4,'linewidth',2);
hold on
%%
Infect5=[];
time=0:0.5:150;
x0=[]; %time初始值（即0的时候）对应的dx=[dE;dI]中E的初始值和I的初始值
S0=300;I0=200;
x0(1)=S0*q;
x0(2)=I0*q;
x0(3)=S0*(1-q);
x0(4)=I0*(1-q);
[t,x]=ode45(@simplify_fun1,time,x0,[]);
for i=1:length(t)
    sum=0;
    for j=2:2:length(x(1,:))
        sum = sum+x(i,j);
    end
    Infect5(i)=sum;
end
h5=plot(t,Infect5,'linewidth',2);
hold on
%%
Infect6=[];
time=0:0.5:150;
x0=[]; %time初始值（即0的时候）对应的dx=[dE;dI]中E的初始值和I的初始值
S0=250;I0=250;
x0(1)=S0*q;
x0(2)=I0*q;
x0(3)=S0*(1-q);
x0(4)=I0*(1-q);
[t,x]=ode45(@simplify_fun1,time,x0,[]);
for i=1:length(t)
    sum=0;
    for j=2:2:length(x(1,:))
        sum = sum+x(i,j);
    end
    Infect6(i)=sum;
end
h6=plot(t,Infect6,'linewidth',2);
hold on
%%
Infect7=[];
time=0:0.5:150;
x0=[]; %time初始值（即0的时候）对应的dx=[dE;dI]中E的初始值和I的初始值
S0=200;I0=300;
x0(1)=S0*q;
x0(2)=I0*q;
x0(3)=S0*(1-q);
x0(4)=I0*(1-q);
[t,x]=ode45(@simplify_fun1,time,x0,[]);
for i=1:length(t)
    sum=0;
    for j=2:2:length(x(1,:))
        sum = sum+x(i,j);
    end
    Infect7(i)=sum;
end
h7=plot(t,Infect7,'linewidth',2);
hold on
%%
Infect8=[];
time=0:0.5:150;
x0=[]; %time初始值（即0的时候）对应的dx=[dE;dI]中E的初始值和I的初始值
S0=150;I0=350;
x0(1)=S0*q;
x0(2)=I0*q;
x0(3)=S0*(1-q);
x0(4)=I0*(1-q);
[t,x]=ode45(@simplify_fun1,time,x0,[]);
for i=1:length(t)
    sum=0;
    for j=2:2:length(x(1,:))
        sum = sum+x(i,j);
    end
    Infect8(i)=sum;
end
h8=plot(t,Infect8,'linewidth',2);
hold on
%%
Infect9=[];
time=0:0.5:150;
x0=[]; %time初始值（即0的时候）对应的dx=[dE;dI]中E的初始值和I的初始值
S0=100;I0=400;
x0(1)=S0*q;
x0(2)=I0*q;
x0(3)=S0*(1-q);
x0(4)=I0*(1-q);
[t,x]=ode45(@simplify_fun1,time,x0,[]);
for i=1:length(t)
    sum=0;
    for j=2:2:length(x(1,:))
        sum = sum+x(i,j);
    end
    Infect9(i)=sum;
end
h9=plot(t,Infect9,'linewidth',2);
hold on
%%
Infect10=[];
time=0:0.5:150;
x0=[]; %time初始值（即0的时候）对应的dx=[dE;dI]中E的初始值和I的初始值
S0=50;I0=450;
x0(1)=S0*q;
x0(2)=I0*q;
x0(3)=S0*(1-q);
x0(4)=I0*(1-q);
[t,x]=ode45(@simplify_fun1,time,x0,[]);
for i=1:length(t)
    sum=0;
    for j=2:2:length(x(1,:))
        sum = sum+x(i,j);
    end
    Infect10(i)=sum;
end
h10=plot(t,Infect10,'linewidth',2);
hold on


%%
xlabel('t');
ylabel('I');
%设置背景透明
set(gca,'color','none');


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