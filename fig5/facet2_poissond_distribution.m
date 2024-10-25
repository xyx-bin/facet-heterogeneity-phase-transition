clc;
clear;
N=2000;
%  %  %  %----------------------------------------泊松分布------------------------------------------------
n=6;   %二阶平均度
degrees_2 = poissrnd(n,N,1);   %服从泊松分布的随机数，n为给定时间间隔内的平均事件数（可以理解为平均度），从n的泊松分布中生成N行1列的随机数
K2=max(degrees_2);                %取degrees中最大的数赋值给K1

%%
% 可视化
figure;
histogram(degrees_2, 'Normalization', 'probability');
k2 = unique(degrees_2);
degree2_Distribution = histcounts(degrees_2,[k2;max(k2)+1],'Normalization','probability');
%---------------------------------------k2----------------------------------------
uniqueDegree_2_ult = 0:max(k2);% 创建一个包含所有可能值的向量
existing_values_2 = ismember(uniqueDegree_2_ult, k2);% 找到向量中存在的值
% 找到缺失的值
missing_values_2 = uniqueDegree_2_ult(~existing_values_2);
missing_probabilities_2 = zeros(1, numel(missing_values_2));
% 在缺失的位置插入 0 元素和相应的概率
for i = 1:numel(missing_values_2)
    index = missing_values_2(i) + 1; % MATLAB索引从1开始，所以需要+1
    degree2_Distribution = [degree2_Distribution(1:index-1) missing_probabilities_2(i) degree2_Distribution(index:end)];
end
average_2=0;
for i=1:length(uniqueDegree_2_ult)
    average_2=average_2+uniqueDegree_2_ult(i)*degree2_Distribution(i);
end
variance_k2=0;
for i=1:length(uniqueDegree_2_ult)
    variance_k2=variance_k2+uniqueDegree_2_ult(i)^2*degree2_Distribution(i);
end
save facet2_poissond_distribution