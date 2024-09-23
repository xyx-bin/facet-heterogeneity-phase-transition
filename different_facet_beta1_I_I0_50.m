clear;
clc;
global beta1 beta2 gamma %定义全局变量（必须，否则会出错）
gamma=0.15;beta2=0.0625;
Beta11=[0:0.0001:0.01];%给定beta1要取的范围

%%
load('facet1=poissrnd_facet2=Zipf_1.mat');
Infect_11=zeros(1,length(Beta11));%每个beta1都会对应一个I，所以需要给I申请一个空间来存放
%通过循环将每个beta1对应的I算出来并放在Infect中
for i=1:length(Beta11)
    beta1=Beta11(i);
    M1=length(uniqueDegree_1_ult);%一维面的长度
    N1=length(uniqueDegree_2_ult);%二维面的长度
    X0=zeros(1,2*M1*N1);
    S0=1950; I0=50;%初始时刻所有易感者和所有感染者的数量
    for k=1:M1 %i代表一维面
        for j=1:N1 %j代表二维面
            X0(2*(k-1)*N1+2*j-1)=S0*Com_num(j,k);
            X0(2*(k-1)*N1+2*j)=I0*Com_num(j,k);
        end
    end
    [t,x]=ode45(@distribution_fun,[0 200],X0,[]);
    sum=0;
    for j=2:2:length(x(1,:))
        sum = sum+x(end,j);
    end
    Infect_11(i)=sum;
    fprintf('beta1=%f时完成！\n',Beta11(i));   
end


%%
load('facet1=poissrnd_facet2=Zipf_2.mat');
Infect_12=zeros(1,length(Beta11));%每个beta1都会对应一个I，所以需要给I申请一个空间来存放
%通过循环将每个beta1对应的I算出来并放在Infect中
for i=1:length(Beta11)
    beta1=Beta11(i);
    M1=length(uniqueDegree_1_ult);%一维面的长度
    N1=length(uniqueDegree_2_ult);%二维面的长度
    X0=zeros(1,2*M1*N1);
    S0=1950; I0=50;%初始时刻所有易感者和所有感染者的数量
    for k=1:M1 %i代表一维面
        for j=1:N1 %j代表二维面
            X0(2*(k-1)*N1+2*j-1)=S0*Com_num(j,k);
            X0(2*(k-1)*N1+2*j)=I0*Com_num(j,k);
        end
    end
    [t,x]=ode45(@distribution_fun2,[0 200],X0,[]);
    sum=0;
    for j=2:2:length(x(1,:))
        sum = sum+x(end,j);
    end
    Infect_12(i)=sum;
    fprintf('beta1=%f时完成！\n',Beta11(i));   
end


%%
load('facet1=poissrnd_facet2=Zipf_3.mat');
Infect_13=zeros(1,length(Beta11));%每个beta1都会对应一个I，所以需要给I申请一个空间来存放
%通过循环将每个beta1对应的I算出来并放在Infect中
for i=1:length(Beta11)
    beta1=Beta11(i);
    M1=length(uniqueDegree_1_ult);%一维面的长度
    N1=length(uniqueDegree_2_ult);%二维面的长度
    X0=zeros(1,2*M1*N1);
    S0=1950; I0=50;%初始时刻所有易感者和所有感染者的数量
    for k=1:M1 %i代表一维面
        for j=1:N1 %j代表二维面
            X0(2*(k-1)*N1+2*j-1)=S0*Com_num(j,k);
            X0(2*(k-1)*N1+2*j)=I0*Com_num(j,k);
        end
    end
    [t,x]=ode45(@distribution_fun3,[0 200],X0,[]);
    sum=0;
    for j=2:2:length(x(1,:))
        sum = sum+x(end,j);
    end
    Infect_13(i)=sum;
    fprintf('beta1=%f时完成！\n',Beta11(i));   
end


%%
load('facet1=poissrnd_facet2=Zipf_4.mat');
Infect_14=zeros(1,length(Beta11));%每个beta1都会对应一个I，所以需要给I申请一个空间来存放
%通过循环将每个beta1对应的I算出来并放在Infect中
for i=1:length(Beta11)
    beta1=Beta11(i);
    M1=length(uniqueDegree_1_ult);%一维面的长度
    N1=length(uniqueDegree_2_ult);%二维面的长度
    X0=zeros(1,2*M1*N1);
    S0=1950; I0=50;%初始时刻所有易感者和所有感染者的数量
    for k=1:M1 %i代表一维面
        for j=1:N1 %j代表二维面
            X0(2*(k-1)*N1+2*j-1)=S0*Com_num(j,k);
            X0(2*(k-1)*N1+2*j)=I0*Com_num(j,k);
        end
    end
    [t,x]=ode45(@distribution_fun4,[0 200],X0,[]);
    sum=0;
    for j=2:2:length(x(1,:))
        sum = sum+x(end,j);
    end
    Infect_14(i)=sum;
    fprintf('beta1=%f时完成！\n',Beta11(i));   
end
h1=plot(Beta11,Infect_11,'Color', '#55B7E6','linewidth',2);
hold on
h2=plot(Beta11,Infect_12,'Color', '#193E8F','linewidth',2);
hold on
h3=plot(Beta11,Infect_13,'Color', '#E53528','linewidth',2);
hold on
h4=plot(Beta11,Infect_14,'Color', '#F09739','linewidth',2);
hold on

%%
legend([h1,h2,h3,h4],'1-facet：power law  2-facet：Poisson','1-facet：power law  2-facet：power law','1-facet：Poisson  2-facet：power law','1-facet：Poisson  2-facet：Poisson','Location','NorthWest');
xlabel('\beta_1');
ylabel('I');

save different_facet1_beta1_I_I0=50
