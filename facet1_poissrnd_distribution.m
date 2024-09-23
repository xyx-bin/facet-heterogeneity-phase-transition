clc;
clear;
N=2000;
%  %  %  %----------------------------------------泊松分布------------------------------------------------
n=20;   %一阶平均度
degrees_1 = poissrnd(n,N,1);   %服从泊松分布的随机数，n为给定时间间隔内的平均事件数（可以理解为平均度），从n的泊松分布中生成N行1列的随机数
K1=max(degrees_1);                %取degrees中最大的数赋值给K1
k1 = unique(degrees_1);
degree1_Distribution = histcounts(degrees_1,[k1;max(k1)+1],'Normalization','probability');
%可视化
figure;
histogram(degrees_1,'Normalization','probability');
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
average_1=0;
for i=1:length(uniqueDegree_1_ult)
    average_1=average_1+uniqueDegree_1_ult(i)*degree1_Distribution(i);
end
variance_k1=0;
for i=1:length(uniqueDegree_1_ult)
    variance_k1=variance_k1+uniqueDegree_1_ult(i)^2*degree1_Distribution(i);
end
save facet1_poissrnd_distribution