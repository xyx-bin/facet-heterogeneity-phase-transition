clc;
clear;
global beta1 beta2 gamma
load('Degree_Distribution.mat');
gamma=0.15;

load('Degree_Distribution.mat');
M=max(uniqueDegree_1);%一阶facet的最大度
N=max(uniqueDegree_2);%二阶facet的最大度
M1=length(uniqueDegree_1);%一维面的长度
N1=length(uniqueDegree_2);%二维面的长度
averagek_1=mean(facet1_degree(:));%一阶facet的平均度
averagek_2=mean(facet2_degree(:));%二阶facet的平均度

%beta1=0.005625;
Beta1 = [0:0.001:0.01];
R=[];
for b = 1:length(Beta1)
    beta1 = Beta1(b);
    X0=zeros(1,2*M1*N1);
    S0=1900; I0=100;%初始时刻所有易感者和所有感染者的密度

    for i=1:M1 %i代表一维面
        for j=1:N1 %j代表二维面
            X0(2*(i-1)*N1+2*j-1)=S0*Com_num(j,i);
            X0(2*(i-1)*N1+2*j)=I0*Com_num(j,i);
        end
    end
    [T1,X1]=ode45(@SIS_fun, [0 300] ,X0);%基本再生数是在无病平衡点处求得的，所以要稳定以后再对N(k1,k2)进行计算

    all = X1(end,:);
    F = [];
    for i=1:M1
        F1 = [];
        for j=1:M1
            N_f=0;
            for l=1:N1
                N_f = N_f + all((i-1)*2*N1+2*l-1) + all((i-1)*2*N1+2*l);
            end  %正确
            N_ff = 0;
            if N_f == 0
                N_ff = 0;
            else
                N_ff = 1/N_f;
            end
            FF=[];
            for k=1:N1
                FF(k) = all((i-1)*2*N1+2*k-1) + all((i-1)*2*N1+2*k);
            end
            F111 = repmat(FF, N1, 1)'; %正确
            F112 = [];
            for s = 1:N1
                N_t = 0;
                N_t = all((i-1)*2*N1+2*s-1) + all((i-1)*2*N1+2*s);
                for t = 1:N1
                    N_n = 0;
                    for k=1:M1
                        N_n=N_n+all((k-1)*2*N1+2*t)+all((k-1)*2*N1+2*t-1);
                    end
                    N_nn = 0;
                    if N_n == 0
                        N_nn = 0;
                    else
                        N_nn = N_n;
                    end
                    if N_nn == 0
                        F112(s,t)=0;
                    else
                        F112(s,t) = (s*N_t*t*sum(Com_num(t,:)))/(N_nn);
                    end                    
                end
            end
            F1 = [F1,beta1*i*(1/averagek_1)*j*sum(Com_num(:,j))*N_ff*F111 + 2*beta1*(1/averagek_2)*F112];
        end
        F = [F;F1];
    end
    Vinv = 1/gamma * diag(ones(M1*N1, 1));
    % FV^-1
    FV = F * Vinv;
    % 计算乘积矩阵的特征值
    eigenvalues = eig(FV);
    % 选择最大特征值
    R0 = max(eigenvalues);
    disp(Beta1(b));
    disp(R0);
    R(b) = R0;
end
h1 = plot(Beta1,R,'-ob','linewidth',1);
hold on;




%%
mu=0.15;
k1=20;
Beta21 = [0:0.001:0.01];
R21=[];
for m=1:length(Beta21)
    R21(m)=(Beta21(m)*k1)/mu;
end
h4=plot(Beta21,R21,'-or','linewidth',1);
hold on;


legend([h1,h2,h3,h4,h5,h6],'<k^2>/<k>=22.87','<k^2>/<k>=27.59','<k^2>/<k>=42.88','<k^2>/<k>=22.87','<k^2>/<k>=27.59','<k^2>/<k>=42.88','\beta_1=0.008','Location','NorthWest');

xlabel('\beta_1');
ylabel('R_0');
set(gca,'color','none');
hold on;

save R0_beta1o