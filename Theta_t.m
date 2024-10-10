clear;
clc;
%%
global beta1 beta2 gamma;

%% 
beta1 = 0.002;% 调0.002、0.008
X0=[0.1,0.1];%初值应该设置为多少才能与fig4a对应
[T1,X1]=ode45(@selfequation1, [0:2:200] ,X0);

%%
beta1 = 0.002;% 调0.002、0.008
X0=[0.4,0.4];%初值应该设置为多少才能与fig4a对应
[T2,X2]=ode45(@selfequation1, [0:2:200] ,X0);

%%
beta1 = 0.008;% 调0.002、0.008
X0=[0.1,0.1];%初值应该设置为多少才能与fig4a对应
[T3,X3]=ode45(@selfequation1, [0:2:200] ,X0);

%%
beta1 = 0.008;% 调0.002、0.008
X0=[0.4,0.4];%初值应该设置为多少才能与fig4a对应
[T4,X4]=ode45(@selfequation1, [0:2:200] ,X0);

%% Theta1
plot(T1,X1(:,1),'-b','linewidth',1);
hold on;
plot(T2,X2(:,1),'-r','linewidth',1);
hold on;
plot(T3,X3(:,1),'-b','linewidth',1);
hold on;
plot(T4,X4(:,1),'-r','linewidth',1);
hold on;
xlabel('t');
ylabel('\Theta_1^*');
set(gca,'color','none');

%%
figure;
plot(T1,X1(:,2),'-b','linewidth',1);
hold on;
plot(T2,X2(:,2),'-r','linewidth',1);
hold on;
plot(T3,X3(:,2),'-b','linewidth',1);
hold on;
plot(T4,X4(:,2),'-r','linewidth',1);
hold on;
xlabel('t');
ylabel('\Theta_2^*');
set(gca,'color','none');


%%
function F=selfequation1(t,x)
%定义参数
global beta1 beta2 gamma;
beta2=0.0625;%会影响（调β_2）β_2不会影响阈值，但是会影响出现双稳的点，所以通过调β_2直到出现双稳，再调β_1
gamma=0.15;%会影响阈值

load('Degree_Distribution.mat');

M=max(uniqueDegree_1_ult);% 1-维面的最大度
N=max(uniqueDegree_2_ult);% 2-维面的最大度
%%
% 1-维面和2-维面的平均度
average_1=0;
for i=1:length(uniqueDegree_1)
    average_1=average_1+uniqueDegree_1(i)*degree1_Distribution(i);
end
average_2=0;
for i=1:length(uniqueDegree_2)
    average_2=average_2+uniqueDegree_2(i)*degree2_Distribution(i);
end

%%
F=zeros(2,1); %给2个方程分配空间
F1 = 0;
for i=1:length(uniqueDegree_1)
    F11 = 0;
    for j=1:length(uniqueDegree_2)
        F14 = Com_num(j,i)*(beta1*uniqueDegree_1(i)*x(1)+2*beta1*uniqueDegree_2(j)*x(2)+beta2*uniqueDegree_2(j)*x(2)^2)/(beta1*uniqueDegree_1(i)*x(1)+2*beta1*uniqueDegree_2(j)*x(2)+beta2*uniqueDegree_2(j)*x(2)^2+gamma);
        if isnan(F14)
            F14 = 0;
        end
        F11 = F11 + F14;
    end
    F12 = 0;
    for k=1:length(uniqueDegree_2)
        F12 = F12 + Com_num(k,i);
    end
    F13 = (uniqueDegree_1(i)*degree1_Distribution(i)) * F11 / F12;
    if isnan(F13)
        F13 = 0;
    end
    F1 = F1 + F13;
end
F(1)=F1/average_1-x(1);

F2=0;
for j=1:length(uniqueDegree_2)
    F21 = 0;
    for i=1:length(uniqueDegree_1)
        F23 = Com_num(j,i)*(beta1*uniqueDegree_1(i)*x(1)+2*beta1*uniqueDegree_2(j)*x(2)+beta2*uniqueDegree_2(j)*x(2)^2)/(beta1*uniqueDegree_1(i)*x(1)+2*beta1*uniqueDegree_2(j)*x(2)+beta2*uniqueDegree_2(j)*x(2)^2+gamma);
        if isnan(F23)
            F23 = 0;
        end
        F21 = F21 + F23;
    end
    F22 = 0;
    for k=1:length(uniqueDegree_1)
        F22 = F22 + Com_num(j,k);
    end
    F24 = (uniqueDegree_2(j)*degree2_Distribution(j)) * F21 / F22;
    if isnan(F24)
        F24 = 0;
    end
    F2 = F2 + F24;
end
F(2)=F2/average_2-x(2);
fprintf('第 %d 时刻的数据计算完成\n', t);
end