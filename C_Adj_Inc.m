%邻接矩阵和关联矩阵的转换
function W = C_Adj_Inc(F,f)
%当f=0时，邻接矩阵转换为关联矩阵，F表示邻接矩阵，W表示关联矩阵
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
%当f=1时，关联矩阵转换为邻接矩阵，F表示关联矩阵，W表示邻接矩阵
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
