function f=SIS(t,x)
%定义参数，全部趋于正平衡点
beta1=0.005625;
beta2=0.0625;%不影响分支
gamma=0.15;

load('Degree_Distribution.mat');

M=max(uniqueDegree_1);%一阶facet的最大度
N=max(uniqueDegree_2);%二阶facet的最大度

%一维面的长度
M1=length(uniqueDegree_1);
%二维面的长度
N1=length(uniqueDegree_2);

averagek_1=mean(facet1_degree(:));%一阶facet的平均度
averagek_2=mean(facet2_degree(:));%二阶facet的平均度

f=zeros(2*M1*N1,1);%给2*M*N个方程分配空间

%%
%计算Θ_1
sum_1=0;
for i=1:M1
    I_s=0;%I(i,1)+I(i,2)+...+I(i,N1)
    for j=1:N1 
        I_s=I_s+x(2*(i-1)*N1+2*j);
    end
    N_s=0;%I(i,1)+I(i,2)+...+I(i,N)+S(i,1)+S(i,2)+...+S(i,N)=N(i,1)+N(i,2)+...+N(i,N)
    for k=1:N1
        N_s=N_s+x(2*(i-1)*N1+2*k)+x(2*(i-1)*N1+2*k-1);   
    end
    sum_1=sum_1+(uniqueDegree_1(i)*degree1_Distribution(i))*(I_s/N_s);
end
theta_1=sum_1/averagek_1;
%计算Θ_2
sum_2=0;
for i=1:N1
    I_n=0; %I(1,i)+I(2,i)+...+I(M,i)
    for j=1:M1
        I_n=I_n+x((j-1)*2*N1+2*i);
    end
    N_n=0;
    for k=1:M1
        N_n=N_n+x((k-1)*2*N1+2*i)+x((k-1)*2*N1+2*i-1);        
    end
    sum_2=sum_2+(uniqueDegree_2(i)*degree2_Distribution(i))*(I_n/N_n);
end
theta_2=sum_2/averagek_2;
%%
fprintf('第%f时刻',t)
for i=1:M1
    for j=1:N1
        f(2*(i-1)*N1+2*j-1) = -beta1*uniqueDegree_1(i)*theta_1*x(2*(i-1)*N1+2*j-1)-2*beta1*uniqueDegree_2(j)*theta_2*x(2*(i-1)*N1+2*j-1)-beta2*uniqueDegree_2(j)*theta_2^2*x(2*(i-1)*N1+2*j-1)+gamma*x(2*(i-1)*N1+2*j);
        f(2*(i-1)*N1+2*j) =  beta1*uniqueDegree_1(i)*theta_1*x(2*(i-1)*N1+2*j-1)+2*beta1*uniqueDegree_2(j)*theta_2*x(2*(i-1)*N1+2*j-1)+beta2*uniqueDegree_2(j)*theta_2^2*x(2*(i-1)*N1+2*j-1)-gamma*x(2*(i-1)*N1+2*j);
    end
end
fprintf('2MN维的数据已生成!\n')



