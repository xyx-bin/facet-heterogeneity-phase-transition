function f=distribution_fun(t,x)
global beta1 beta2 gamma
beta2=0.0625;%不影响分支
gamma=0.15;

load('facet1=Zipf_facet2=poissrnd_1.mat');
%一维面的长度
M1=length(uniqueDegree_1_ult);
%二维面的长度
N1=length(uniqueDegree_2_ult);

averagek_1=mean(degrees_1(:));%一阶facet的平均度
averagek_2=mean(degrees_2(:));%二阶facet的平均度

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
    if N_s == 0
        sum_divide_1 = 0;
    else 
        sum_divide_1 = I_s/N_s;
    end
    sum_1=sum_1+(uniqueDegree_1_ult(i)*degree1_Distribution(i))*(sum_divide_1);
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
    if N_s == 0
        sum_divide_2 = 0;
    else 
        sum_divide_2 = I_s/N_s;
    end
    sum_2=sum_2+(uniqueDegree_2_ult(i)*degree2_Distribution(i))*(sum_divide_2);
end
theta_2=sum_2/averagek_2;
%%
%fprintf('第%f时刻',t)
for i=1:M1
    for j=1:N1
        f(2*(i-1)*N1+2*j-1) = -beta1*uniqueDegree_1_ult(i)*theta_1*x(2*(i-1)*N1+2*j-1)-2*beta1*uniqueDegree_2_ult(j)*theta_2*x(2*(i-1)*N1+2*j-1)-beta2*uniqueDegree_2_ult(j)*theta_2^2*x(2*(i-1)*N1+2*j-1)+gamma*x(2*(i-1)*N1+2*j);
        f(2*(i-1)*N1+2*j) =  beta1*uniqueDegree_1_ult(i)*theta_1*x(2*(i-1)*N1+2*j-1)+2*beta1*uniqueDegree_2_ult(j)*theta_2*x(2*(i-1)*N1+2*j-1)+beta2*uniqueDegree_2_ult(j)*theta_2^2*x(2*(i-1)*N1+2*j-1)-gamma*x(2*(i-1)*N1+2*j);
    end
end
%fprintf('2MN维的数据已生成!\n')



