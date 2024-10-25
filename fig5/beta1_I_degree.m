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
plot(Beta11',Infect1,'-r');
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
plot(Beta12',Infect2,'-r');
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

plot(Beta13',Infect3,'-r');
hold on;

%%
clc;
clear;
load('Theta1_date_2287.mat');
load('Theta2_date_2287.mat');
load('facet1=Zipf_facet2=poissrnd_2.mat');
global beta1;

beta2=0.0625;%会影响（调β_2）β_2不会影响阈值，但是会影响出现双稳的点，所以通过调β_2直到出现双稳，再调β_1
gamma=0.15;%会影响阈值

Infect4=[];
Beta1=Beta31;%给定beta1要取的范围
for k=1:length(Beta1)
    beta1=Beta1(k);
    I=0;
    for i=1:length(uniqueDegree_1_ult)
        for j=1:length(uniqueDegree_2_ult)
            I = I + (Com_num(j,i)*(beta1*uniqueDegree_1_ult(i)*Theta31(k)+2*beta1*uniqueDegree_2_ult(j)*Theta41(k)+beta2*uniqueDegree_2_ult(j)*Theta41(k)^2))/(beta1*uniqueDegree_1_ult(i)*Theta31(k)+2*beta1*uniqueDegree_2_ult(j)*Theta41(k)+beta2*uniqueDegree_2_ult(j)*Theta41(k)^2+gamma);
        end
    end
    Infect4(k) = I;
end
plot(Beta31',Infect4,'-b');
hold on;
%%
Infect5=[];
Beta1=Beta32;%给定beta1要取的范围
for k=1:length(Beta1)
    beta1=Beta1(k);
    I=0;
    for i=1:length(uniqueDegree_1_ult)
        for j=1:length(uniqueDegree_2_ult)
            I = I + (Com_num(j,i)*(beta1*uniqueDegree_1_ult(i)*Theta32(k)+2*beta1*uniqueDegree_2_ult(j)*Theta42(k)+beta2*uniqueDegree_2_ult(j)*Theta42(k)^2))/(beta1*uniqueDegree_1_ult(i)*Theta32(k)+2*beta1*uniqueDegree_2_ult(j)*Theta42(k)+beta2*uniqueDegree_2_ult(j)*Theta42(k)^2+gamma);
        end
    end
    Infect5(k) = I;
end
plot(Beta32',Infect5,'-b');
hold on;

%%
Infect6=[];
Beta1=Beta33;%给定beta1要取的范围
for k=1:length(Beta1)
    beta1=Beta1(k);
    I=0;
    for i=1:length(uniqueDegree_1_ult)
        for j=1:length(uniqueDegree_2_ult)
            I = I + (Com_num(j,i)*(beta1*uniqueDegree_1_ult(i)*Theta33(k)+2*beta1*uniqueDegree_2_ult(j)*Theta43(k)+beta2*uniqueDegree_2_ult(j)*Theta43(k)^2))/(beta1*uniqueDegree_1_ult(i)*Theta33(k)+2*beta1*uniqueDegree_2_ult(j)*Theta43(k)+beta2*uniqueDegree_2_ult(j)*Theta43(k)^2+gamma);
        end
    end
    Infect6(k) = I;
end

plot(Beta33',Infect6,'-b');
hold on;

%%
clc;
clear;
load('Theta1_date_2759.mat');
load('Theta2_date_2759.mat');
load('facet1=Zipf_facet2=poissrnd_3.mat');
global beta1;

beta2=0.0625;%会影响（调β_2）β_2不会影响阈值，但是会影响出现双稳的点，所以通过调β_2直到出现双稳，再调β_1
gamma=0.15;%会影响阈值

Infect7=[];
Beta1=Beta51;%给定beta1要取的范围
for k=1:length(Beta1)
    beta1=Beta1(k);
    I=0;
    for i=1:length(uniqueDegree_1_ult)
        for j=1:length(uniqueDegree_2_ult)
            I = I + (Com_num(j,i)*(beta1*uniqueDegree_1_ult(i)*Theta51(k)+2*beta1*uniqueDegree_2_ult(j)*Theta61(k)+beta2*uniqueDegree_2_ult(j)*Theta61(k)^2))/(beta1*uniqueDegree_1_ult(i)*Theta51(k)+2*beta1*uniqueDegree_2_ult(j)*Theta61(k)+beta2*uniqueDegree_2_ult(j)*Theta61(k)^2+gamma);
        end
    end
    Infect7(k) = I;
end
plot(Beta51',Infect7,'-y');
hold on;
%%
Infect8=[];
Beta1=Beta52;%给定beta1要取的范围
for k=1:length(Beta1)
    beta1=Beta1(k);
    I=0;
    for i=1:length(uniqueDegree_1_ult)
        for j=1:length(uniqueDegree_2_ult)
            I = I + (Com_num(j,i)*(beta1*uniqueDegree_1_ult(i)*Theta52(k)+2*beta1*uniqueDegree_2_ult(j)*Theta62(k)+beta2*uniqueDegree_2_ult(j)*Theta62(k)^2))/(beta1*uniqueDegree_1_ult(i)*Theta52(k)+2*beta1*uniqueDegree_2_ult(j)*Theta62(k)+beta2*uniqueDegree_2_ult(j)*Theta62(k)^2+gamma);
        end
    end
    Infect8(k) = I;
end
plot(Beta52',Infect8,'-y');
hold on;

%%
Infect9=[];
Beta1=Beta53;%给定beta1要取的范围
for k=1:length(Beta1)
    beta1=Beta1(k);
    I=0;
    for i=1:length(uniqueDegree_1_ult)
        for j=1:length(uniqueDegree_2_ult)
            I = I + (Com_num(j,i)*(beta1*uniqueDegree_1_ult(i)*Theta53(k)+2*beta1*uniqueDegree_2_ult(j)*Theta63(k)+beta2*uniqueDegree_2_ult(j)*Theta63(k)^2))/(beta1*uniqueDegree_1_ult(i)*Theta53(k)+2*beta1*uniqueDegree_2_ult(j)*Theta63(k)+beta2*uniqueDegree_2_ult(j)*Theta63(k)^2+gamma);
        end
    end
    Infect9(k) = I;
end

plot(Beta53',Infect9,'-y');
hold on;

%%
clc;
clear;
load('Theta1_date_4288.mat');
load('Theta2_date_4288.mat');
load('facet1=Zipf_facet2=poissrnd_4.mat');
global beta1;

beta2=0.0625;%会影响（调β_2）β_2不会影响阈值，但是会影响出现双稳的点，所以通过调β_2直到出现双稳，再调β_1
gamma=0.15;%会影响阈值

Infect10=[];
Beta1=Beta71;%给定beta1要取的范围
for k=1:length(Beta1)
    beta1=Beta1(k);
    I=0;
    for i=1:length(uniqueDegree_1_ult)
        for j=1:length(uniqueDegree_2_ult)
            I = I + (Com_num(j,i)*(beta1*uniqueDegree_1_ult(i)*Theta71(k)+2*beta1*uniqueDegree_2_ult(j)*Theta81(k)+beta2*uniqueDegree_2_ult(j)*Theta81(k)^2))/(beta1*uniqueDegree_1_ult(i)*Theta71(k)+2*beta1*uniqueDegree_2_ult(j)*Theta81(k)+beta2*uniqueDegree_2_ult(j)*Theta81(k)^2+gamma);
        end
    end
    Infect10(k) = I;
end
plot(Beta71',Infect10,'-g');
hold on;
%%
Infect11=[];
Beta1=Beta72;%给定beta1要取的范围
for k=1:length(Beta1)
    beta1=Beta1(k);
    I=0;
    for i=1:length(uniqueDegree_1_ult)
        for j=1:length(uniqueDegree_2_ult)
            I = I + (Com_num(j,i)*(beta1*uniqueDegree_1_ult(i)*Theta72(k)+2*beta1*uniqueDegree_2_ult(j)*Theta82(k)+beta2*uniqueDegree_2_ult(j)*Theta82(k)^2))/(beta1*uniqueDegree_1_ult(i)*Theta72(k)+2*beta1*uniqueDegree_2_ult(j)*Theta82(k)+beta2*uniqueDegree_2_ult(j)*Theta82(k)^2+gamma);
        end
    end
    Infect11(k) = I;
end
plot(Beta72',Infect11,'-g');
hold on;

%%
Infect12=[];
Beta1=Beta73;%给定beta1要取的范围
for k=1:length(Beta1)
    beta1=Beta1(k);
    I=0;
    for i=1:length(uniqueDegree_1_ult)
        for j=1:length(uniqueDegree_2_ult)
            I = I + (Com_num(j,i)*(beta1*uniqueDegree_1_ult(i)*Theta73(k)+2*beta1*uniqueDegree_2_ult(j)*Theta83(k)+beta2*uniqueDegree_2_ult(j)*Theta83(k)^2))/(beta1*uniqueDegree_1_ult(i)*Theta73(k)+2*beta1*uniqueDegree_2_ult(j)*Theta83(k)+beta2*uniqueDegree_2_ult(j)*Theta83(k)^2+gamma);
        end
    end
    Infect12(k) = I;
end

plot(Beta73',Infect12,'-g');
hold on;


xlabel('\beta_1');
ylabel('I');