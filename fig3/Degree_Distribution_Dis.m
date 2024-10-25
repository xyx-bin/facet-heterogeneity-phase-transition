clear; clc
N=2000;
k1 = 20;   %一阶平均度
k2 = 6;%二阶平均度
%----------------------生成一阶、二阶单纯形的邻接矩阵（根据一阶A1、二阶平均度生成的邻接矩阵A2）------------------------------------------
p1 = zeros(length(k2),1);   %存放两个节点之间连接的概率的25*1的矩阵
p2 = zeros(length(k2),1);

for m = 1:1:length(k2)
    k2 = k2(m);
    p1 = (k1-2*k2)/((N-1)-2*k2);%任意两个节点之间连接的概率
    p2 = 2 * k2 / ((N-1)*(N-2));%任意三个节点构成三角形的概率
    
    p1(m) = p1;
    p2(m) = p2;
    
    %一阶单纯形。存放根据p1生成的矩阵（根据一阶平均度生成的一阶单纯形的邻接矩阵A1）
    A1 = zeros(N,N);
    for i = 1:N
        for j =i+1:N
            if rand<p1%产生(0,1)均匀分布中随机得到的随机数
                A1(i,j)=1;%上三角矩阵(1,2),(1,3),(1,4)......(2,3),(2,4)......(N-1,N)
                A1(j,i)=1;%下三角矩阵(2,1),(3,1),(4,1)......(3,2),(4,2)......(N,N-1)
            end
        end
    end
    
    %二阶单纯形。根据二阶平均度生成的邻接矩阵A2
    count_triad = 0;
    triad = [];
    for i =1:N
        for j = i+1:N
            for k = j+1:N
                if rand<p2
                    A1(i,j) = 1;
                    A1(j,i) = 1;
                    A1(i,k) = 1;
                    A1(k,i) = 1;
                    A1(k,j) = 1;
                    A1(j,k) = 1;
                    triad = [triad;i,j,k];%记录组成三角形的节点(i,j,k)的位置，并且是以i节点从小到大的顺序排列的(没有重复)
                    count_triad = count_triad+1;%三角形的数量
                end
            end
        end
    end
    
    %计算网络上二阶单纯形的邻接矩阵
    incidence_2 = zeros(N,size(triad,1)); %  incidence matrix: nodes belonging to which triad (order 2/ 2-simplices edge)
    %size(triad,1)返回triad矩阵的行数(三角形的个数)；size(triad,2)返回triad矩阵的列数；
    for i =1:N
        [e,~] = find(triad==i);
        incidence_2(i,e) = 1;
    end%第i个节点连接的三角形的个数(只记了一次)？？？？？？？？？
    A2 = incidence_2 * incidence_2';
    A2 = A2 - diag(diag(A2)); %二阶邻接矩阵，确实是diag了两次，一次的话对角线上有元素？？？？？？？？？？
end
C=A1+A2;%A1是一阶邻接矩阵，A2是二阶邻接矩阵；C是单纯复形网络的邻接矩阵


%%
%----------------------------------------------二阶facet的度分布----------------------------------------------------
%每个节点连接的三角形（二阶facet）的个数，即为关联矩阵的行和
facet2_degree = sum(incidence_2~=0,2);%每行非零元素的个数
uniqueDegree_2 = unique(facet2_degree);%找出唯一的元素（每个节点连接三角形的数量）
%计算度分布，Normalization设置为probability表示一概率形式计数
degree2_Distribution = histcounts(facet2_degree,[uniqueDegree_2;max(uniqueDegree_2)+1],'Normalization','probability');
%绘制度分布
bar(uniqueDegree_2,degree2_Distribution);
xlabel('facet2 degree');
ylabel('degree2 Distribution');
hold on;

%%
missing_values_2 = [];
%找出uniqueDegree_2中缺失的值
for i = 2:length(uniqueDegree_2)
    if uniqueDegree_2(i) - uniqueDegree_2(i-1) > 1
        missing_values_2 = [missing_values_2, uniqueDegree_2(i-1)+1:uniqueDegree_2(i)-1];
    end
end
if (uniqueDegree_2(1)==0)
    uniqueDegree_2_ult=uniqueDegree_2(2:end);
    degree2_Distribution_ult = degree2_Distribution(2:end);
end
%如果缺失，则在degree2_Distribution缺失的位置处添加0
if (~isempty(missing_values_2))
    for i=1:length(missing_values_2)
        uniqueDegree_2_ult=[uniqueDegree_2_ult(1:missing_values_2(i)-1); missing_values_2(i); uniqueDegree_2_ult(missing_values_2(i):end)];
        degree2_Distribution_ult = [degree2_Distribution_ult(1:missing_values_2(i)-1), 0, degree2_Distribution_ult(missing_values_2(i):end)];
    end
end


%%
%----------------------------------------------一阶facet的度分布----------------------------------------------------
%查找连接节点i的三角形中共用一条边的个数
repeat_edges_N = zeros(N,1);
for i=1:N
    repeat_num=[0,0];
    inc_i_triad=[];
    [line_i,column_i]=find(triad==i);
    inc_i_triad=[inc_i_triad;triad(line_i,:)];
    for j=2:size(inc_i_triad,1)
        for k=1:j-1
            num = intersect(inc_i_triad(j,:),inc_i_triad(k,:));
            if (length(num)==2)
                if(k==1)%每添加一个三角形第一次比较的时候，但凡有公共边的就会计数+1；第二次以后比较有公共边则不+（否则会重复计数）；
                    repeat_edges_N(i) = repeat_edges_N(i)+1;
                    repeat_num = [repeat_num;num];
                else if(ismember(num, repeat_num, 'rows')~=1)
                        repeat_edges_N(i) = repeat_edges_N(i)+1;
                        repeat_num = [repeat_num;num];
                    end
                end
            end
        end
    end
end

% 单纯复形网络中节点i的度为第i行非零元素的个数，而不是第i行元素相加
incidence_1 = C_Adj_Inc(A1,0);%计算一阶邻接矩阵的关联矩阵
facet1_degree = sum(incidence_1~=0,2);%每行非零元素的个数
%每个节点连接的单边的数量（不包括三角形中的边）
for i=1:size(facet1_degree,1)
    facet1_degree(i) =  facet1_degree(i) + repeat_edges_N(i) - facet2_degree(i)*2;
end
uniqueDegree_1 = unique(facet1_degree);%找出唯一的元素
degree1_Distribution = histcounts(facet1_degree,[uniqueDegree_1;max(uniqueDegree_1)+1],'Normalization','probability');
figure;
bar(uniqueDegree_1,degree1_Distribution);
xlabel('facet1 degree');
ylabel('degree1 Distribution');
%%
missing_values_1 = [];
%找出uniqueDegree_2中缺失的值
for i = 2:length(uniqueDegree_1)
    if uniqueDegree_1(i) - uniqueDegree_1(i-1) > 1
        missing_values_1 = [missing_values_1, uniqueDegree_1(i-1)+1:uniqueDegree_1(i)-1];
    end
end
if (uniqueDegree_1(1)==0)
    uniqueDegree_1_ult=uniqueDegree_1(2:end);
    degree1_Distribution_ult = degree1_Distribution(2:end);
end
%如果缺失，则在degree2_Distribution缺失的位置处添加0
if (~isempty(missing_values_1))
    for i=1:length(missing_values_1)
        uniqueDegree_1_ult=[uniqueDegree_1_ult(1:missing_values_1(i)-1); missing_values_1(i); uniqueDegree_1_ult(missing_values_1(i):end)];
        degree1_Distribution_ult = [degree1_Distribution_ult(1:missing_values_1(i)-1), 0, degree1_Distribution_ult(missing_values_1(i):end)];
    end
end

%%
%计算一维面和二维面的联合概率分布
Com_num = zeros(length(uniqueDegree_2),length(uniqueDegree_1));%网络的联合概率分布
for i=0:length(uniqueDegree_1)-1
    res = find(facet1_degree==i);
    temp = zeros(1,length(res));
    for k=1:length(res)
        temp(k)=facet2_degree(res(k));
    end
    for j=0:length(uniqueDegree_2)-1
        num=0;
        if(find(temp==j))
            num = num + histcounts(temp,[j-0.5,j+0.5]);%计算数组中某个数字出现的次数
        end
        Com_num(j+1,i+1)=num/N;
    end
end

save Degree_Distribution



