
clc;
clear;
load('Degree_Distribution.mat');

% 计算每个节点的度
degrees = sum(A1, 2);

% 对节点的度进行排序，并统计相同度数的节点数
[sorted_degrees, ~, idx] = unique(degrees);
node_counts = zeros(1, length(sorted_degrees));
for i = 1:length(sorted_degrees)
    node_counts(i) = sum(degrees(idx) == sorted_degrees(i));
end

% 将节点数放在一个矩阵中
node_count_matrix = node_counts;

% 计算每个度数出现的概率
probabilities = node_counts / sum(node_counts);

% 将概率放在另一个矩阵中
probability_matrix = probabilities;

% 显示结果
disp('从小到大连接相同数量的节点的度:');
disp(sorted_degrees);
disp('对应的概率矩阵:');
disp(probability_matrix);

average_degree=0;
for i=1:length(sorted_degrees)    
    average_degree = average_degree + sorted_degrees(i)*probability_matrix(i);
end
disp('网络的平均度');
disp(average_degree);

second_order_degree = 0;
for i=1:length(sorted_degrees)    
    second_order_degree = second_order_degree + sorted_degrees(i)^2*probability_matrix(i);
end
disp('网络的二阶矩');
disp(second_order_degree);

heterogeneity=0;
heterogeneity = second_order_degree/average_degree;
disp('网络的异质性');
disp(heterogeneity);