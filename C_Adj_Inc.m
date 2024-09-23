%�ڽӾ���͹��������ת��
function W = C_Adj_Inc(F,f)
%��f=0ʱ���ڽӾ���ת��Ϊ��������F��ʾ�ڽӾ���W��ʾ��������
if f == 0
    m = sum(sum(F))/2;
    n = size(F,1);
    W = zeros(n,m);
    k = 1;
    for i = 1 : n
        for j = i : n
            if(F(i,j) ~= 0)
                W(i,k) = 1;
                W(j,k) = 1;
                k = k + 1;
            end
        end
    end
end
%��f=1ʱ����������ת��Ϊ�ڽӾ���F��ʾ��������W��ʾ�ڽӾ���
if f == 1
    m = size(F,2);
    n = size(F,1);
    W = zeros(n,n);
    for i = 1 : m
        a = find(F(:,i) ~= 0);
        W(a(1),a(2)) = 1;
        W(a(2),a(1)) = 1;
    end
end
