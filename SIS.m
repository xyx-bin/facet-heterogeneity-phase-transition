function f=SIS(t,x)
%���������ȫ��������ƽ���
beta1=0.005625;
beta2=0.0625;%��Ӱ���֧
gamma=0.15;

load('Degree_Distribution.mat');

M=max(uniqueDegree_1);%һ��facet������
N=max(uniqueDegree_2);%����facet������

%һά��ĳ���
M1=length(uniqueDegree_1);
%��ά��ĳ���
N1=length(uniqueDegree_2);

averagek_1=mean(facet1_degree(:));%һ��facet��ƽ����
averagek_2=mean(facet2_degree(:));%����facet��ƽ����

f=zeros(2*M1*N1,1);%��2*M*N�����̷���ռ�

%%
%���㦨_1
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
%���㦨_2
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
fprintf('��%fʱ��',t)
for i=1:M1
    for j=1:N1
        f(2*(i-1)*N1+2*j-1) = -beta1*uniqueDegree_1(i)*theta_1*x(2*(i-1)*N1+2*j-1)-2*beta1*uniqueDegree_2(j)*theta_2*x(2*(i-1)*N1+2*j-1)-beta2*uniqueDegree_2(j)*theta_2^2*x(2*(i-1)*N1+2*j-1)+gamma*x(2*(i-1)*N1+2*j);
        f(2*(i-1)*N1+2*j) =  beta1*uniqueDegree_1(i)*theta_1*x(2*(i-1)*N1+2*j-1)+2*beta1*uniqueDegree_2(j)*theta_2*x(2*(i-1)*N1+2*j-1)+beta2*uniqueDegree_2(j)*theta_2^2*x(2*(i-1)*N1+2*j-1)-gamma*x(2*(i-1)*N1+2*j);
    end
end
fprintf('2MNά������������!\n')



