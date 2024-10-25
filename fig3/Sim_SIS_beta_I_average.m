clc;
clear;
load('data_beta1_N_2000_I_50_500_U.mat');
h1=plot(Beta11,Infect_11,'-b','linewidth',2);
hold on
h2=plot(Beta1,Infect,'-.r','linewidth',1);
hold on

load('dataChaosSimSIS_I_50_1.mat');
load('dataChaosSimSIS_I_50_2.mat');
load('dataChaosSimSIS_I_50_3.mat');
load('dataChaosSimSIS_I_50_4.mat');
load('dataChaosSimSIS_I_50_5.mat');
load('dataChaosSimSIS_I_50_6.mat');
load('dataChaosSimSIS_I_50_7.mat');
load('dataChaosSimSIS_I_50_8.mat');
load('dataChaosSimSIS_I_50_9.mat');
load('dataChaosSimSIS_I_50_10.mat');
load('dataChaosSimSIS_I_50_11.mat');
load('dataChaosSimSIS_I_50_12.mat');
load('dataChaosSimSIS_I_50_13.mat');
load('dataChaosSimSIS_I_50_14.mat');
load('dataChaosSimSIS_I_50_15.mat');
load('dataChaosSimSIS_I_50_16.mat');
load('dataChaosSimSIS_I_50_17.mat');
load('dataChaosSimSIS_I_50_18.mat');
load('dataChaosSimSIS_I_50_19.mat');
load('dataChaosSimSIS_I_50_20.mat');
Infect_average_50=[];
Infect_average_50=mean([Infect_I_50_1;Infect_I_50_2;Infect_I_50_3;Infect_I_50_4;Infect_I_50_5;Infect_I_50_6;Infect_I_50_7;Infect_I_50_8;Infect_I_50_9;Infect_I_50_10;Infect_I_50_11;Infect_I_50_12;Infect_I_50_13;Infect_I_50_14;Infect_I_50_15;Infect_I_50_16;Infect_I_50_17;Infect_I_50_18;Infect_I_50_19;Infect_I_50_20]);
h3=plot(Beta1,Infect_average_50,'b*','LineWidth',3);

load('dataChaosSimSIS_I_500_1.mat');
load('dataChaosSimSIS_I_500_2.mat');
load('dataChaosSimSIS_I_500_3.mat');
load('dataChaosSimSIS_I_500_4.mat');
load('dataChaosSimSIS_I_500_5.mat');
load('dataChaosSimSIS_I_500_6.mat');
load('dataChaosSimSIS_I_500_7.mat');
load('dataChaosSimSIS_I_500_8.mat');
load('dataChaosSimSIS_I_500_9.mat');
load('dataChaosSimSIS_I_500_10.mat');
load('dataChaosSimSIS_I_500_11.mat');
load('dataChaosSimSIS_I_500_12.mat');
load('dataChaosSimSIS_I_500_13.mat');
load('dataChaosSimSIS_I_500_14.mat');
load('dataChaosSimSIS_I_500_15.mat');
load('dataChaosSimSIS_I_500_16.mat');
load('dataChaosSimSIS_I_500_17.mat');
load('dataChaosSimSIS_I_500_18.mat');
load('dataChaosSimSIS_I_500_19.mat');
load('dataChaosSimSIS_I_500_20.mat');
load('dataChaosSimSIS_I_500_21.mat');
load('dataChaosSimSIS_I_500_22.mat');
load('dataChaosSimSIS_I_500_23.mat');
load('dataChaosSimSIS_I_500_24.mat');
load('dataChaosSimSIS_I_500_25.mat');
load('dataChaosSimSIS_I_500_26.mat');
load('dataChaosSimSIS_I_500_27.mat');
load('dataChaosSimSIS_I_500_28.mat');
load('dataChaosSimSIS_I_500_29.mat');
load('dataChaosSimSIS_I_500_30.mat');
load('dataChaosSimSIS_I_500_31.mat');
load('dataChaosSimSIS_I_500_32.mat');
load('dataChaosSimSIS_I_500_33.mat');
load('dataChaosSimSIS_I_500_34.mat');
load('dataChaosSimSIS_I_500_35.mat');
load('dataChaosSimSIS_I_500_36.mat');
load('dataChaosSimSIS_I_500_37.mat');
load('dataChaosSimSIS_I_500_38.mat');
load('dataChaosSimSIS_I_500_39.mat');
load('dataChaosSimSIS_I_500_40.mat');
Infect_average_500=[];
Infect_average_500=mean([Infect_I_500_1;Infect_I_500_2;Infect_I_500_3;Infect_I_500_4;Infect_I_500_5;Infect_I_500_6;Infect_I_500_7;Infect_I_500_8;Infect_I_500_9;Infect_I_500_10; ...
    Infect_I_500_11;Infect_I_500_12;Infect_I_500_13;Infect_I_500_14;Infect_I_500_15;Infect_I_500_16;Infect_I_500_17;Infect_I_500_18;Infect_I_500_19;Infect_I_500_20; ...
    Infect_I_500_21;Infect_I_500_22;Infect_I_500_23;Infect_I_500_24;Infect_I_500_25;Infect_I_500_26;Infect_I_500_27;Infect_I_500_28;Infect_I_500_29;Infect_I_500_30; ...
    Infect_I_500_31;Infect_I_500_32;Infect_I_500_33;Infect_I_500_34;Infect_I_500_35;Infect_I_500_36;Infect_I_500_37;Infect_I_500_38;Infect_I_500_39;Infect_I_500_40]);
h4=plot(Beta1,Infect_average_500,'r*','LineWidth',3);%Ëæ»úÄ£Äâ
legend([h1,h2,h3,h4],'I=50','I=500','MC I=50','MC I=500','Location','NorthWest');
xlabel('\beta_1');
ylabel('I');
%ÉèÖÃ±³¾°Í¸Ã÷
set(gca,'color','none');


