clc;
clear;
global beta1 beta2 k1 k2 gamma q theta2 N average_2 % 定义全局变量
gamma=0.6;beta2=0.55;N=500;average_2=3; % 这儿的k1,k2是指一维面度为k1，二维面度为k2
% k1=4;
% q=0.1;
% k2=round((average_2-k1*q)/(1-q)); % round四舍五入取整

%%
Beta10=[];
xlab=[];
k1 = 1;%%%%%%%调
q = 0.2;%%%%%%%调
k2=round((average_2-k1*q)/(1-q)); % round四舍五入取整
xlab(1)=(k1^2*q+k2^2*(1-q))/(k1*q+k2*(1-q)); % k^2/k    %%%%%%调
Infect1=[];
Beta1=[0:0.001:0.3];
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
end
plot(Beta1,Infect1,'linewidth',2);
hold on;
Infect2=[];
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
    Infect2(i)=x(end,2)+x(end,4);
end
% plot(Beta1,Infect2,'linewidth',2);
% hold on
if Infect1(1) > 100
    beta10 = 0;
else
    for i=2:length(Infect1)
        if Infect1(i)-Infect1(i-1) > 3
            beta10 = Beta1(i-1);
            break;
        end
    end
end
if Infect2(1) > 100
    beta11 = 0;
else
    for i=2:length(Infect2)
        if Infect2(i)-Infect2(i-1) > 3
            beta11 = Beta1(i-1);
            break;
        end
    end
end
Beta10(1)=min(beta10,beta11);%%%%%%%调

%%
k1 = 2;
q = 0.2;
k2=round((average_2-k1*q)/(1-q)); % round四舍五入取整
xlab(2)=(k1^2*q+k2^2*(1-q))/(k1*q+k2*(1-q)); % k^2/k
Infect1=[];
Beta1=[0:0.001:0.3];
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
end
% plot(Beta1,Infect1,'linewidth',2);
% hold on;
Infect2=[];
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
    Infect2(i)=x(end,2)+x(end,4);
end
% plot(Beta1,Infect2,'linewidth',2);
% hold on
if Infect1(1) > 100
    beta10 = 0;
else
    for i=2:length(Infect1)
        if Infect1(i)-Infect1(i-1) > 3
            beta10 = Beta1(i-1);
            break;
        end
    end
end
if Infect2(1) > 100
    beta11 = 0;
else
    for i=2:length(Infect2)
        if Infect2(i)-Infect2(i-1) > 3
            beta11 = Beta1(i-1);
            break;
        end
    end
end
Beta10(2)=min(beta10,beta11);
%%
k1 = 1;%%%%%%%调
q = 0.1;%%%%%%%调
k2=round((average_2-k1*q)/(1-q)); % round四舍五入取整
xlab(3)=(k1^2*q+k2^2*(1-q))/(k1*q+k2*(1-q)); % k^2/k    %%%%%%调
Infect1=[];
Beta1=[0:0.001:0.3];
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
end
% plot(Beta1,Infect1,'linewidth',2);
% hold on;
Infect2=[];
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
    Infect2(i)=x(end,2)+x(end,4);
end
% plot(Beta1,Infect2,'linewidth',2);
% hold on
if Infect1(1) > 100
    beta10 = 0;
else
    for i=2:length(Infect1)
        if Infect1(i)-Infect1(i-1) > 3
            beta10 = Beta1(i-1);
            break;
        end
    end
end
if Infect2(1) > 100
    beta11 = 0;
else
    for i=2:length(Infect2)
        if Infect2(i)-Infect2(i-1) > 3
            beta11 = Beta1(i-1);
            break;
        end
    end
end
Beta10(3)=min(beta10,beta11);%%%%%%%调

%%
k1 = 2;
q = 0.1;
k2=round((average_2-k1*q)/(1-q)); % round四舍五入取整
xlab(4)=(k1^2*q+k2^2*(1-q))/(k1*q+k2*(1-q)); % k^2/k
Infect1=[];
Beta1=[0:0.001:0.3];
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
end
% plot(Beta1,Infect1,'linewidth',2);
% hold on;
Infect2=[];
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
    Infect2(i)=x(end,2)+x(end,4);
end
% plot(Beta1,Infect2,'linewidth',2);
% hold on
if Infect1(1) > 100
    beta10 = 0;
else
    for i=2:length(Infect1)
        if Infect1(i)-Infect1(i-1) > 3
            beta10 = Beta1(i-1);
            break;
        end
    end
end
if Infect2(1) > 100
    beta11 = 0;
else
    for i=2:length(Infect2)
        if Infect2(i)-Infect2(i-1) > 3
            beta11 = Beta1(i-1);
            break;
        end
    end
end
Beta10(4)=min(beta10,beta11);

%%
k1 = 3;
q = 0.1;
k2=round((average_2-k1*q)/(1-q)); % round四舍五入取整
xlab(5)=(k1^2*q+k2^2*(1-q))/(k1*q+k2*(1-q)); % k^2/k
Infect1=[];
Beta1=[0:0.001:0.3];
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
end
% plot(Beta1,Infect1,'linewidth',2);
% hold on;
Infect2=[];
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
    Infect2(i)=x(end,2)+x(end,4);
end
% plot(Beta1,Infect2,'linewidth',2);
% hold on
if Infect1(1) > 100
    beta10 = 0;
else
    for i=2:length(Infect1)
        if Infect1(i)-Infect1(i-1) > 3
            beta10 = Beta1(i-1);
            break;
        end
    end
end
if Infect2(1) > 100
    beta11 = 0;
else
    for i=2:length(Infect2)
        if Infect2(i)-Infect2(i-1) > 3
            beta11 = Beta1(i-1);
            break;
        end
    end
end
Beta10(5)=min(beta10,beta11);

%%
k1 =4;
q = 0.1;
k2=round((average_2-k1*q)/(1-q)); % round四舍五入取整
xlab(6)=(k1^2*q+k2^2*(1-q))/(k1*q+k2*(1-q)); % k^2/k
Infect1=[];
Beta1=[0:0.001:0.3];
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
end
% plot(Beta1,Infect1,'linewidth',2);
% hold on;
Infect2=[];
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
    Infect2(i)=x(end,2)+x(end,4);
end
% plot(Beta1,Infect2,'linewidth',2);
% hold on
if Infect1(1) > 100
    beta10 = 0;
else
    for i=2:length(Infect1)
        if Infect1(i)-Infect1(i-1) > 3
            beta10 = Beta1(i-1);
            break;
        end
    end
end
if Infect2(1) > 100
    beta11 = 0;
else
    for i=2:length(Infect2)
        if Infect2(i)-Infect2(i-1) > 3
            beta11 = Beta1(i-1);
            break;
        end
    end
end
Beta10(6)=min(beta10,beta11);

%%
k1 = 5;
q = 0.1;
k2=round((average_2-k1*q)/(1-q)); % round四舍五入取整
xlab(7)=(k1^2*q+k2^2*(1-q))/(k1*q+k2*(1-q)); % k^2/k
Infect1=[];
Beta1=[0:0.001:0.3];
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
end
% plot(Beta1,Infect1,'linewidth',2);
% hold on;
Infect2=[];
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
    Infect2(i)=x(end,2)+x(end,4);
end
% plot(Beta1,Infect2,'linewidth',2);
% hold on
if Infect1(1) > 100
    beta10 = 0;
else
    for i=2:length(Infect1)
        if Infect1(i)-Infect1(i-1) > 3
            beta10 = Beta1(i-1);
            break;
        end
    end
end
if Infect2(1) > 100
    beta11 = 0;
else
    for i=2:length(Infect2)
        if Infect2(i)-Infect2(i-1) > 3
            beta11 = Beta1(i-1);
            break;
        end
    end
end
Beta10(7)=min(beta10,beta11);

%%
k1 = 6;
q = 0.1;
k2=round((average_2-k1*q)/(1-q)); % round四舍五入取整
xlab(8)=(k1^2*q+k2^2*(1-q))/(k1*q+k2*(1-q)); % k^2/k
Infect1=[];
Beta1=[0:0.001:0.3];
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
end
% plot(Beta1,Infect1,'linewidth',2);
% hold on;
Infect2=[];
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
    Infect2(i)=x(end,2)+x(end,4);
end
% plot(Beta1,Infect2,'linewidth',2);
% hold on
if Infect1(1) > 100
    beta10 = 0;
else
    for i=2:length(Infect1)
        if Infect1(i)-Infect1(i-1) > 3
            beta10 = Beta1(i-1);
            break;
        end
    end
end
if Infect2(1) > 100
    beta11 = 0;
else
    for i=2:length(Infect2)
        if Infect2(i)-Infect2(i-1) > 3
            beta11 = Beta1(i-1);
            break;
        end
    end
end
Beta10(8)=min(beta10,beta11);

%%
k1 = 7;
q = 0.1;
k2=round((average_2-k1*q)/(1-q)); % round四舍五入取整
xlab(9)=(k1^2*q+k2^2*(1-q))/(k1*q+k2*(1-q)); % k^2/k
Infect1=[];
Beta1=[0:0.001:0.3];
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
end
% plot(Beta1,Infect1,'linewidth',2);
% hold on;
Infect2=[];
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
    Infect2(i)=x(end,2)+x(end,4);
end
% plot(Beta1,Infect2,'linewidth',2);
% hold on
if Infect1(1) > 100
    beta10 = 0;
else
    for i=2:length(Infect1)
        if Infect1(i)-Infect1(i-1) > 3
            beta10 = Beta1(i-1);
            break;
        end
    end
end
if Infect2(1) > 100
    beta11 = 0;
else
    for i=2:length(Infect2)
        if Infect2(i)-Infect2(i-1) > 3
            beta11 = Beta1(i-1);
            break;
        end
    end
end
Beta10(9)=min(beta10,beta11);



%%
save q_beta10


plot(xlab,Beta10,'linewidth',2);
hold on
xlim([average_2 9]);
ylim([0 0.3]);


%%
%function [y1,...,yN] = myfun(x1,...,xM) .名为 myfun 的函数，
%该函数接受输入 x1,...,xM 并返回输出 y1,...,yN
function dx=simplify_fun1(t,x)
%修改变量为全局变量
global beta1 beta2 k1 k2 gamma q N theta2 average_2
dx=zeros(4,1);
theta2=(k1*x(2)+k2*x(4))/(N*average_2);
dx(1)=-2*beta1*k1*x(1)*theta2-beta2*k1*x(1)*theta2^2+gamma*x(2);
dx(2)=2*beta1*k1*x(1)*theta2+beta2*k1*x(1)*theta2^2-gamma*x(2);
dx(3)=-2*beta1*k2*x(3)*theta2-beta2*k2*x(3)*theta2^2+gamma*x(4);
dx(4)=2*beta1*k2*x(3)*theta2+beta2*k2*x(3)*theta2^2-gamma*x(4);
end