clear; clc
N=2000;
k1 = 20;   %һ��ƽ����
k2 = 6;%����ƽ����
%----------------------����һ�ס����׵����ε��ڽӾ��󣨸���һ��A1������ƽ�������ɵ��ڽӾ���A2��------------------------------------------
p1 = zeros(length(k2),1);   %��������ڵ�֮�����ӵĸ��ʵ�25*1�ľ���
p2 = zeros(length(k2),1);

for m = 1:1:length(k2)
    k2 = k2(m);
    p1 = (k1-2*k2)/((N-1)-2*k2);%���������ڵ�֮�����ӵĸ���
    p2 = 2 * k2 / ((N-1)*(N-2));%���������ڵ㹹�������εĸ���
    
    p1(m) = p1;
    p2(m) = p2;
    
    %һ�׵����Ρ���Ÿ���p1���ɵľ��󣨸���һ��ƽ�������ɵ�һ�׵����ε��ڽӾ���A1��
    A1 = zeros(N,N);
    for i = 1:N
        for j =i+1:N
            if rand<p1%����(0,1)���ȷֲ�������õ��������
                A1(i,j)=1;%�����Ǿ���(1,2),(1,3),(1,4)......(2,3),(2,4)......(N-1,N)
                A1(j,i)=1;%�����Ǿ���(2,1),(3,1),(4,1)......(3,2),(4,2)......(N,N-1)
            end
        end
    end
    
    %���׵����Ρ����ݶ���ƽ�������ɵ��ڽӾ���A2
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
                    triad = [triad;i,j,k];%��¼��������εĽڵ�(i,j,k)��λ�ã���������i�ڵ��С�����˳�����е�(û���ظ�)
                    count_triad = count_triad+1;%�����ε�����
                end
            end
        end
    end
    
    %���������϶��׵����ε��ڽӾ���
    incidence_2 = zeros(N,size(triad,1)); %  incidence matrix: nodes belonging to which triad (order 2/ 2-simplices edge)
    %size(triad,1)����triad���������(�����εĸ���)��size(triad,2)����triad�����������
    for i =1:N
        [e,~] = find(triad==i);
        incidence_2(i,e) = 1;
    end%��i���ڵ����ӵ������εĸ���(ֻ����һ��)������������������
    A2 = incidence_2 * incidence_2';
    A2 = A2 - diag(diag(A2)); %�����ڽӾ���ȷʵ��diag�����Σ�һ�εĻ��Խ�������Ԫ�أ�������������������
end
C=A1+A2;%A1��һ���ڽӾ���A2�Ƕ����ڽӾ���C�ǵ�������������ڽӾ���


%%
%----------------------------------------------����facet�Ķȷֲ�----------------------------------------------------
%ÿ���ڵ����ӵ������Σ�����facet���ĸ�������Ϊ����������к�
facet2_degree = sum(incidence_2~=0,2);%ÿ�з���Ԫ�صĸ���
uniqueDegree_2 = unique(facet2_degree);%�ҳ�Ψһ��Ԫ�أ�ÿ���ڵ����������ε�������
%����ȷֲ���Normalization����Ϊprobability��ʾһ������ʽ����
degree2_Distribution = histcounts(facet2_degree,[uniqueDegree_2;max(uniqueDegree_2)+1],'Normalization','probability');
%���ƶȷֲ�
bar(uniqueDegree_2,degree2_Distribution);
xlabel('facet2 degree');
ylabel('degree2 Distribution');
hold on;

%%
missing_values_2 = [];
%�ҳ�uniqueDegree_2��ȱʧ��ֵ
for i = 2:length(uniqueDegree_2)
    if uniqueDegree_2(i) - uniqueDegree_2(i-1) > 1
        missing_values_2 = [missing_values_2, uniqueDegree_2(i-1)+1:uniqueDegree_2(i)-1];
    end
end
if (uniqueDegree_2(1)==0)
    uniqueDegree_2_ult=uniqueDegree_2(2:end);
    degree2_Distribution_ult = degree2_Distribution(2:end);
end
%���ȱʧ������degree2_Distributionȱʧ��λ�ô����0
if (~isempty(missing_values_2))
    for i=1:length(missing_values_2)
        uniqueDegree_2_ult=[uniqueDegree_2_ult(1:missing_values_2(i)-1); missing_values_2(i); uniqueDegree_2_ult(missing_values_2(i):end)];
        degree2_Distribution_ult = [degree2_Distribution_ult(1:missing_values_2(i)-1), 0, degree2_Distribution_ult(missing_values_2(i):end)];
    end
end


%%
%----------------------------------------------һ��facet�Ķȷֲ�----------------------------------------------------
%�������ӽڵ�i���������й���һ���ߵĸ���
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
                if(k==1)%ÿ���һ�������ε�һ�αȽϵ�ʱ�򣬵����й����ߵľͻ����+1���ڶ����Ժ�Ƚ��й�������+��������ظ���������
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

% �������������нڵ�i�Ķ�Ϊ��i�з���Ԫ�صĸ����������ǵ�i��Ԫ�����
incidence_1 = C_Adj_Inc(A1,0);%����һ���ڽӾ���Ĺ�������
facet1_degree = sum(incidence_1~=0,2);%ÿ�з���Ԫ�صĸ���
%ÿ���ڵ����ӵĵ��ߵ��������������������еıߣ�
for i=1:size(facet1_degree,1)
    facet1_degree(i) =  facet1_degree(i) + repeat_edges_N(i) - facet2_degree(i)*2;
end
uniqueDegree_1 = unique(facet1_degree);%�ҳ�Ψһ��Ԫ��
degree1_Distribution = histcounts(facet1_degree,[uniqueDegree_1;max(uniqueDegree_1)+1],'Normalization','probability');
figure;
bar(uniqueDegree_1,degree1_Distribution);
xlabel('facet1 degree');
ylabel('degree1 Distribution');
%%
missing_values_1 = [];
%�ҳ�uniqueDegree_2��ȱʧ��ֵ
for i = 2:length(uniqueDegree_1)
    if uniqueDegree_1(i) - uniqueDegree_1(i-1) > 1
        missing_values_1 = [missing_values_1, uniqueDegree_1(i-1)+1:uniqueDegree_1(i)-1];
    end
end
if (uniqueDegree_1(1)==0)
    uniqueDegree_1_ult=uniqueDegree_1(2:end);
    degree1_Distribution_ult = degree1_Distribution(2:end);
end
%���ȱʧ������degree2_Distributionȱʧ��λ�ô����0
if (~isempty(missing_values_1))
    for i=1:length(missing_values_1)
        uniqueDegree_1_ult=[uniqueDegree_1_ult(1:missing_values_1(i)-1); missing_values_1(i); uniqueDegree_1_ult(missing_values_1(i):end)];
        degree1_Distribution_ult = [degree1_Distribution_ult(1:missing_values_1(i)-1), 0, degree1_Distribution_ult(missing_values_1(i):end)];
    end
end

%%
%����һά��Ͷ�ά������ϸ��ʷֲ�
Com_num = zeros(length(uniqueDegree_2),length(uniqueDegree_1));%��������ϸ��ʷֲ�
for i=0:length(uniqueDegree_1)-1
    res = find(facet1_degree==i);
    temp = zeros(1,length(res));
    for k=1:length(res)
        temp(k)=facet2_degree(res(k));
    end
    for j=0:length(uniqueDegree_2)-1
        num=0;
        if(find(temp==j))
            num = num + histcounts(temp,[j-0.5,j+0.5]);%����������ĳ�����ֳ��ֵĴ���
        end
        Com_num(j+1,i+1)=num/N;
    end
end

save Degree_Distribution



