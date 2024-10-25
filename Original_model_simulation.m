%ctrl+Rһ��ע�Ͷ���
%ctrl+Tȡ��ע��
%ctrl+C��ֹ
%���΢�ַ������һ�ַ�����4���弶Runge-Kutta�㷨
%[t,x]=ode45(Fun,tspan,x0,options,pars)
%t����������XΪn��n��������ÿһ�о�������һ��״̬������t�ı仯ֵ
%funΪҪ���ķ����飬��ʽ����Ϊһ��΢�ַ����飻@���������ұߵĲ��֣�������������
%tspanΪҪ����t�����䣬x0Ϊ��ʼֵ��options��һЩѡ��
clear;
clc;
global beta1 beta2 gamma %����ȫ�ֱ��������룬��������
gamma=0.15;beta2=0.0625;
load('Degree_Distribution.mat');

%%
Beta11=[0:0.001:0.03];%����beta1Ҫȡ�ķ�Χ
Infect_11=zeros(1,length(Beta11));%ÿ��beta1�����Ӧһ��I��������Ҫ��I����һ���ռ������
%ͨ��ѭ����ÿ��beta1��Ӧ��I�����������Infect��
for i=1:length(Beta11)
    beta1=Beta11(i);
    M1=length(uniqueDegree_1);%һά��ĳ���
    N1=length(uniqueDegree_2);%��ά��ĳ���
    X0=zeros(1,2*M1*N1);
    S0=1950; I0=50;%��ʼʱ�������׸��ߺ����и�Ⱦ�ߵ�����
    for k=1:M1 %i����һά��
        for j=1:N1 %j�����ά��
            X0(2*(k-1)*N1+2*j-1)=S0*Com_num(j,k);
            X0(2*(k-1)*N1+2*j)=I0*Com_num(j,k);
        end
    end
    [t,x]=ode45(@SIS_fun,[0 225],X0,[]);
    sum=0;
    for j=2:2:length(x(1,:))
        sum = sum+x(end,j);
    end
    Infect_11(i)=sum;
    fprintf('beta1=%fʱ��ɣ�\n',Beta11(i));   
end
h1=plot(Beta11,Infect_11,'-b','linewidth',2);
hold on
%%
Beta1=[0:0.001:0.03];%����beta1Ҫȡ�ķ�Χ
Infect=zeros(1,length(Beta1));%ÿ��beta1�����Ӧһ��I��������Ҫ��I����һ���ռ������
%ͨ��ѭ����ÿ��beta1��Ӧ��I�����������Infect��
for i=1:length(Beta1)
    beta1=Beta1(i);
    M1=length(uniqueDegree_1);%һά��ĳ���
    N1=length(uniqueDegree_2);%��ά��ĳ���
    X0=zeros(1,2*M1*N1);
    S0=1500; I0=500;%��ʼʱ�������׸��ߺ����и�Ⱦ�ߵ��ܶ�
    for k=1:M1 %i����һά��
        for j=1:N1 %j�����ά��
            X0(2*(k-1)*N1+2*j-1)=S0*Com_num(j,k);
            X0(2*(k-1)*N1+2*j)=I0*Com_num(j,k);
        end
    end
    [t,x]=ode45(@SIS_fun,[0 225],X0,[]);
    sum=0;
    for j=2:2:length(x(1,:))
        sum = sum+x(end,j);
    end
    Infect(i)=sum;
    fprintf('beta1=%fʱ��ɣ�\n',Beta1(i));
end
h2=plot(Beta1,Infect,'-.r','linewidth',1);
hold on

%%
%ͼ����"\phi"�������ϣ����ĸ��_��������±꣬\infty�����������
legend([h1,h2],'I=50','I=500','Location','NorthWest');
xlabel('\beta_1');
ylabel('I');
%���ñ���͸��
set(gca,'color','none');
save data_beta1_N_2000_I_50_500_U





