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

save Theta1_2759
%%
function F=selfequation1(x)
%定义参数
global beta1;
% beta1=0.02;%会影响
beta2=0.0625;%会影响（调β_2）β_2不会影响阈值，但是会影响出现双稳的点，所以通过调β_2直到出现双稳，再调β_1
gamma=0.15;%会影响阈值

load('facet1=Zipf_facet2=poissrnd_3.mat');

M=max(uniqueDegree_1_ult);% 1-维面的最大度
N=max(uniqueDegree_2_ult);% 2-维面的最大度


%%
% 1-维面和2-维面的平均度
average_1=0;
for i=1:length(uniqueDegree_1_ult)
    average_1=average_1+uniqueDegree_1_ult(i)*degree1_Distribution(i);
end
average_2=0;
for i=1:length(uniqueDegree_2_ult)
    average_2=average_2+uniqueDegree_2_ult(i)*degree2_Distribution(i);
end

%%
F=zeros(2,1); %给2个方程分配空间
F1 = 0;
for i=1:length(uniqueDegree_1_ult)
    F11 = 0;
    for j=1:length(uniqueDegree_2_ult)
        F14 = Com_num(j,i)*(beta1*uniqueDegree_1_ult(i)*x(1)+2*beta1*uniqueDegree_2_ult(j)*x(2)+beta2*uniqueDegree_2_ult(j)*x(2)^2)/(beta1*uniqueDegree_1_ult(i)*x(1)+2*beta1*uniqueDegree_2_ult(j)*x(2)+beta2*uniqueDegree_2_ult(j)*x(2)^2+gamma);
        if isnan(F14)
            F14 = 0;
        end
        F11 = F11 + F14;
    end
    F12 = 0;
    for k=1:length(uniqueDegree_2_ult)
        F12 = F12 + Com_num(k,i);
    end
    F13 = (uniqueDegree_1_ult(i)*degree1_Distribution(i)) * F11 / F12;
    if isnan(F13)
        F13 = 0;
    end
    F1 = F1 + F13;
end
F(1)=F1/average_1-x(1);

F2=0;
for j=1:length(uniqueDegree_2_ult)
    F21 = 0;
    for i=1:length(uniqueDegree_1_ult)
        F23 = Com_num(j,i)*(beta1*uniqueDegree_1_ult(i)*x(1)+2*beta1*uniqueDegree_2_ult(j)*x(2)+beta2*uniqueDegree_2_ult(j)*x(2)^2)/(beta1*uniqueDegree_1_ult(i)*x(1)+2*beta1*uniqueDegree_2_ult(j)*x(2)+beta2*uniqueDegree_2_ult(j)*x(2)^2+gamma);
        if isnan(F23)
            F23 = 0;
        end
        F21 = F21 + F23;
    end
    F22 = 0;
    for k=1:length(uniqueDegree_1_ult)
        F22 = F22 + Com_num(j,k);
    end
    F24 = (uniqueDegree_2_ult(j)*degree2_Distribution(j)) * F21 / F22;
    if isnan(F24)
        F24 = 0;
    end
    F2 = F2 + F24;
end
F(2)=F2/average_2-x(2);

end
