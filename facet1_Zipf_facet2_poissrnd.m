clc;
clear;
N=2000;
% % % %----------------------------幂率分布------没有孤立节点-------生成遵循幂律分布的度序列-------------------------

alpha=3;
K=32;%给定最大度
for i=15:K
    p(i)=i^(-alpha);%幂率分布
end
C=1/sum(p);%标准化
Nk=round(C*p*N);%不同度的节点的数量
%语法形式只有1种：Y = round(X)，这里的X可以是数，向量，矩阵，输出对应;用于舍入到最接近的整数(四舍五入)。
LNk=length(find(Nk~=0));%找出不同度的节点数量不为零的度的个数（长度）
degrees1=[];%定义一个空矩阵
for i=1:length(Nk)
    if Nk(i)~=0
        degrees1=[degrees1,ones(1,Nk(i))*i];%ones(1,Nk(i))生成1行Nk(i)列的全是1的矩阵。
    end
end
degrees1=[degrees1,zeros(1,N-length(degrees1))];%缺失的点数为孤立点
%y = randsample([50:100],20) 返回一个向量，从整数 50 到 100 中无放回随机均匀抽取的 20 个值。
degrees_1=randsample(degrees1,length(degrees1));%生成的幂率度序列（从degrees1中不放回的抽取500个值放入degree中，相当于对degrees1中的数随机打乱）
K1=max(degrees_1);

%%
% 可视化
figure;
histogram(degrees_1, 'Normalization', 'probability');

% %  %  %  %----------------------------------------泊松分布------------------------------------------------
% n=6;   %二阶平均度
% degrees_2 = poissrnd(n,N,1);   %服从泊松分布的随机数，n为给定时间间隔内的平均事件数（可以理解为平均度），从n的泊松分布中生成N行1列的随机数
% K2=max(degrees_2);                %取degrees中最大的数赋值给K1
% 
% %%
% %可视化
% figure;
% histogram(degrees_2, 'Normalization', 'probability');

% % % %----------------------------------------联合概率分布----------------------------------------------
k1 = unique(degrees_1);
% k2 = unique(degrees_2);

%一维面的度分布
degree1_Distribution = histcounts(degrees_1',[k1';max(k1')+1],'Normalization','probability');
%二维面的度分布
% degree2_Distribution = histcounts(degrees_2,[k2;max(k2)+1],'Normalization','probability');

%%
% % % %---------------------------将k1与k2中的缺少的元素用0补全-----------------------------------
%---------------------------------------k1-----------------------------------------
uniqueDegree_1_ult = 0:max(k1);% 创建一个包含所有可能值的向量
existing_values_1 = ismember(uniqueDegree_1_ult, k1);% 找到向量中存在的值
% 找到缺失的值
missing_values_1 = uniqueDegree_1_ult(~existing_values_1);
missing_probabilities_1 = zeros(1, numel(missing_values_1));
% 在缺失的位置插入 0 元素和相应的概率
for i = 1:numel(missing_values_1)
    index = missing_values_1(i) + 1; % MATLAB索引从1开始，所以需要+1
    degree1_Distribution = [degree1_Distribution(1:index-1) missing_probabilities_1(i) degree1_Distribution(index:end)];
end


% %---------------------------------------k2----------------------------------------
% uniqueDegree_2_ult = 0:max(k2);% 创建一个包含所有可能值的向量
% existing_values_2 = ismember(uniqueDegree_2_ult, k2);% 找到向量中存在的值
% % 找到缺失的值
% missing_values_2 = uniqueDegree_2_ult(~existing_values_2);
% missing_probabilities_2 = zeros(1, numel(missing_values_2));
% % 在缺失的位置插入 0 元素和相应的概率
% for i = 1:numel(missing_values_2)
%     index = missing_values_2(i) + 1; % MATLAB索引从1开始，所以需要+1
%     degree2_Distribution = [degree2_Distribution(1:index-1) missing_probabilities_2(i) degree2_Distribution(index:end)];
% end

load('facet2_poissrnd_distribution.mat');
%%
Com_num = zeros(length(uniqueDegree_2_ult),length(uniqueDegree_1_ult));%网络的联合概率分布
for i = 1:N
    j = i;
    Com_num(degrees_2(i)+1,degrees_1(i)+1) = Com_num(degrees_2(i)+1,degrees_1(j)+1) + 1/N;
end

%%
% ---------------------计算真正的平均度----------------------------------------
average_1=0;
for i=1:length(uniqueDegree_1_ult)
    average_1=average_1+uniqueDegree_1_ult(i)*degree1_Distribution(i);
end
% average_2=0;
% for i=1:length(uniqueDegree_2_ult)
%     average_2=average_2+uniqueDegree_2_ult(i)*degree2_Distribution(i);
% end
variance_k1=0;
for i=1:length(uniqueDegree_1_ult)
    variance_k1=variance_k1+uniqueDegree_1_ult(i)^2*degree1_Distribution(i);
end
% variance_k2=0;
% for i=1:length(uniqueDegree_2_ult)
%     variance_k2=variance_k2+uniqueDegree_2_ult(i)^2*degree2_Distribution(i);
% end
xlab=variance_k1/average_1;

save facet1=Zipf_facet2=poissrnd