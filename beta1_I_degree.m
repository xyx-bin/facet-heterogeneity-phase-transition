clc;
clear;

load('Theta1_degree_date.mat');
load('Theta2_degree_date.mat');
load('Degree_Distribution.mat');
global beta1;

beta2=0.0625;%会影响（调β_2）β_2不会影响阈值，但是会影响出现双稳的点，所以通过调β_2直到出现双稳，再调β_1
gamma=0.15;%会影响阈值

Infect1=[];
Beta1=Beta11;%给定beta1要取的范围
for k=1:length(Beta1)
    beta1=Beta1(k);
    I=0;
    for i=1:length(uniqueDegree_1)
        for j=1:length(uniqueDegree_2)
            I = I + (Com_num(j,i)*(beta1*uniqueDegree_1(i)*Theta11(k)+2*beta1*uniqueDegree_2(j)*Theta21(k)+beta2*uniqueDegree_2(j)*Theta21(k)^2))/(beta1*uniqueDegree_1(i)*Theta11(k)+2*beta1*uniqueDegree_2(j)*Theta21(k)+beta2*uniqueDegree_2(j)*Theta21(k)^2+gamma);
        end
    end
    Infect1(k) = I;
end
plot(Beta11',Infect1);
hold on;
%%
Infect2=[];
Beta1=Beta12;%给定beta1要取的范围
for k=1:length(Beta1)
    beta1=Beta1(k);
    I=0;
    for i=1:length(uniqueDegree_1)
        for j=1:length(uniqueDegree_2)
            I = I + (Com_num(j,i)*(beta1*uniqueDegree_1(i)*Theta12(k)+2*beta1*uniqueDegree_2(j)*Theta22(k)+beta2*uniqueDegree_2(j)*Theta22(k)^2))/(beta1*uniqueDegree_1(i)*Theta12(k)+2*beta1*uniqueDegree_2(j)*Theta22(k)+beta2*uniqueDegree_2(j)*Theta22(k)^2+gamma);
        end
    end
    Infect2(k) = I;
end
plot(Beta12',Infect2);
hold on;

%%
Infect3=[];
Beta1=Beta13;%给定beta1要取的范围
for k=1:length(Beta1)
    beta1=Beta1(k);
    I=0;
    for i=1:length(uniqueDegree_1)
        for j=1:length(uniqueDegree_2)
            I = I + (Com_num(j,i)*(beta1*uniqueDegree_1(i)*Theta13(k)+2*beta1*uniqueDegree_2(j)*Theta23(k)+beta2*uniqueDegree_2(j)*Theta23(k)^2))/(beta1*uniqueDegree_1(i)*Theta13(k)+2*beta1*uniqueDegree_2(j)*Theta23(k)+beta2*uniqueDegree_2(j)*Theta23(k)^2+gamma);
        end
    end
    Infect3(k) = I;
end

plot(Beta13',Infect3);

xlabel('\beta_1');
ylabel('I');