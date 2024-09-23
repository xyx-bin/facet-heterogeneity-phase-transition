clear;
clc;
global beta1

%%
Beta1=[0:0.00025:0.03];%给定beta1要取的范围
Infect_1_1=zeros(1,length(Beta1));%每个beta1都会对应一个I，所以需要给I申请一个空间来存放
%通过循环将每个beta1对应的I算出来并放在Infect中
for k=1:length(Beta1)
    beta1=Beta1(k);
    for i=0:0.25:1
        x0 = [0;i]; %theta1 和 theta2 的初始猜测值
        x= fsolve(@selfequation1, x0);
        if x(1)>=0
            scatter(beta1, x(1));
            Infect_1_1=[Infect_1_1,x(1)];
            hold on
        end
    end 
    fprintf('beta1=%f 的点已生成!\n',beta1);
end
%将Infect_1_1中的元素从小到大排序

xlabel('\beta_1');
ylabel('\theta_1');
hold on

save Theta1_degree
%%
function F=selfequation1(x)
%定义参数
global beta1;
% beta1=0.02;%会影响
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
% 1-维面度分布和2-维面度分布的二阶矩
second_moment_1=0;
for i=1:length(uniqueDegree_1)
    second_moment_1=second_moment_1+uniqueDegree_1(i)^2*degree1_Distribution(i);
end
second_moment_2=0;
for j=1:length(uniqueDegree_2)
    second_moment_2=second_moment_2+uniqueDegree_2(j)^2*degree2_Distribution(j);
end

%阈值的计算，网络固定，gamma固定，阈值就固定
%threshold_value = (gamma*averagek_1*averagek_2)/(averagek_2*second_moment_1+2*averagek_1*second_moment_2);

%%
F=zeros(2,1); %给2个方程分配空间
F1 = 0;
for i=1:length(uniqueDegree_1)
    F11 = 0;
    for j=1:length(uniqueDegree_2)
        F11 = F11 + Com_num(j,i)*(beta1*uniqueDegree_1(i)*x(1)+2*beta1*uniqueDegree_2(j)*x(2)+beta2*uniqueDegree_2(j)*x(2)^2)/(beta1*uniqueDegree_1(i)*x(1)+2*beta1*uniqueDegree_2(j)*x(2)+beta2*uniqueDegree_2(j)*x(2)^2+gamma);
    end
    if isnan(F11)
        F11 = 0;
    end
    F12 = 0;
    for k=1:length(uniqueDegree_2)
        F12 = F12 + Com_num(k,i);
    end
    F1 = F1 + (uniqueDegree_1(i)*degree1_Distribution(i)) * F11 / F12;
    if isnan(F1)
        F1 = 0;
    end
end
F(1)=F1/average_1-x(1);

F2=0;
for j=1:length(uniqueDegree_2)
    F21 = 0;
    for i=1:length(uniqueDegree_1)
        F21 = F21 + Com_num(j,i)*(beta1*uniqueDegree_1(i)*x(1)+2*beta1*uniqueDegree_2(j)*x(2)+beta2*uniqueDegree_2(j)*x(2)^2)/(beta1*uniqueDegree_1(i)*x(1)+2*beta1*uniqueDegree_2(j)*x(2)+beta2*uniqueDegree_2(j)*x(2)^2+gamma);
    end
    if isnan(F21)
        F21 = 0;
    end
    F22 = 0;
    for k=1:length(uniqueDegree_1)
        F22 = F22 + Com_num(j,k);
    end
    F2 = F2 + (uniqueDegree_2(j)*degree2_Distribution(j)) * F21 / F22;
    if isnan(F22)
        F22 = 0;
    end
end
F(2)=F2/average_2-x(2);

end
