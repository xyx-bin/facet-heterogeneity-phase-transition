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

    %一阶单纯形。存放根据p1生成的矩阵（根据一阶平均度生成的邻接矩阵A1）
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
    e2 = zeros(N,size(triad,1)); %  incidence matrix: nodes belonging to which triad (order 2/ 2-simplices edge)
    %size(triad,1)返回triad矩阵的行数(三角形的个数)；size(triad,2)返回triad矩阵的列数；
    for i =1:N
        [x,~] = find(triad==i);
        e2(i,x) = 1;
    end%第i个节点连接的三角形的个数(只记了一次)？？？？？？？？？
    A2 = e2 * e2';
    A2 = A2 - diag(diag(A2)); %二阶邻接矩阵，确实是diag了两次，一次的话对角线上有元素？？？？？？？？？？
end

C=A1+A2;%A1是一阶邻接矩阵，A2是二阶邻接矩阵；C是单纯复形网络的邻接矩阵

%===============================SIS疾病传播================================

%--------------------------------参数初值----------------------------------
Beta1=[0:0.001:0.03];
Infect_I_50_1=[];
for m=1:length(Beta1)
    beta=Beta1(m);

    %beta=0.02;                   %beta=tau
    beta2=0.0625;
    gamma=0.15;
    I0=50;
    degrees =full(sum(A1));%把稀疏矩阵转换为全矩阵;sum(A1)是对A1的列求和（每个节点的度）
    N=length(degrees);
    time=0;
    nodes=randsample(N,I0);%返回从整数1到N中无放回随机均匀抽取的I0个值
    susp=ones(1,N);%1行n列的全1的矩阵，所有的易感感染节点
    susp(nodes)=0;

    %找出各种类型的三角形
    i_num=0;%记录三角形中感染节点的个数
    sss_triad=[];
    iss_triad=[];
    iis_triad=[];
    iii_triad=[];
    for i=1:count_triad
        ri_triad=triad(i,:);
        for j=1:length(ri_triad)
            for k=1:length(nodes)
                if (ri_triad(j)==nodes(k))
                    i_num=i_num+1; %如果节点相同，则说明此节点是感染节点，记1，若还有相同的，累加
                end
            end
        end
        if(i_num==0)
            sss_triad=[sss_triad;ri_triad];%SSS
        end
        if(i_num==1)
            iss_triad=[iss_triad;ri_triad];%ISS
        end
        if (i_num==2)
            iis_triad=[iis_triad;ri_triad];%将IIS类型的三角形放入矩阵中
        end
        if (i_num==3)
            iii_triad=[iii_triad;ri_triad];%III
        end
        i_num=0;
    end


    %-----------------寻找SIS模型的初值[SI],[II],[I]---------------------
    s_ind=find(susp==1);%寻找易感节点的地址
    i_ind=find(susp==0);%寻找染病节点的地址
    SS=full(sum(sum(A1(s_ind,s_ind))));%A1(s_ind,s_ind)是两个位置之间有连边
    SI=full(sum(sum(A1(s_ind,i_ind))));
    II=full(sum(sum(A1(i_ind,i_ind))));%II对关系的数量
    ee=[time,SS];%初始时刻对应的SS的数量
    cc=[time,SI];
    gg=[time,II];
    mm=[time,size(iis_triad,1)];%初始时刻对应的IIS三角形的数量
    I=length(i_ind);
    xx=[time, I];

    %---------------------------SIS疾病传播----------------------------------------------------------------
    while time<200 && (~isempty(nodes) || ~isempty(iis_triad) || ~isempty(iii_triad) || ~isempty(iss_triad) || ~isempty(sss_triad))
        %----------------------------疾病传播-------------------------------------------------
        if (isempty(nodes)==1)
            Tinf=Inf;
        else
            for i=1:length(nodes)%nodes是感染节点
                if(degrees(nodes(i))==0)%因为泊松分布中有些度为0，如果第一次抽到这样的节点，它不会传染
                    Tinf(i)=Inf;%无穷大，-Inf是无穷小
                else
                    Tinf(i)=exprnd(1/(degrees(nodes(i))*beta))+time;%感染或者恢复的时间间隔都服从指数分布，exprnd(mu) 从均值为 mu 的指数分布中生成一个随机数
                end
            end
        end

        %按三角形IIS传播的时间序列地址（与上面无关，单独发生）
        Tinf2=[];
        if (size(iis_triad,1)==0)
            Tinf2=Inf;
        else
            for i=1:size(iis_triad,1)
                Tinf2(i) = exprnd(1/(1*beta2))+time;
            end
        end
        %exprnd(mu,sz1,...,szN) 从指数分布中生成一个随机数数组，其中 sz1,...,szN 表示每个维度的大小。
        if (isempty(nodes)==1)
            Trec=Inf;
        else
            Trec=exprnd(1/gamma,length(nodes),1)+time;%恢复时间序列地址
        end
        %找出最小的是因为用时时间最短的先感染和先恢复
        [next_inf , ind_inf]=min(Tinf);%将最小值赋值给next_inf,将最小值对应的i赋值给ind_inf
        [next_rec , ind_rec]=min(Trec);
        [next_inf2 , ind_inf2]=min(Tinf2);%
        if (length(nodes)==0)
            break;
        else
            node_inf = nodes(ind_inf);
            node_rec = nodes(ind_rec);
        end
        if (size(iis_triad,1)==0)
            next_inf2=Inf;%如果没有可以感染的IIS三角形，就取不出时间最短，会出错，所以需要终止
        else
            triad_inf = iis_triad((ind_inf2),:);%时间最短的三角形
        end

        %--------====有一条宗旨：从A添加到B，A中对应的项要删去=====----------------------------------
        if (next_inf < next_rec && next_inf < next_inf2)
            time = next_inf;
            neighbors = find(A1(node_inf,:));%找出传染节点的邻居(在这个邻居当中可能是三角形的邻居，这样的话就会令三角形变成III，但是并未在iis_triad 中删除)
            ind=randsample(neighbors,1);%返回从neighbors中无放回随机均匀抽取的1个值（ind可能为S，也可能为I）
            ind_in_triad=[];%存放新感染节点在iis_triad中所在行的三角形
            [line_number,column_number]=find(iis_triad==ind);
            ind_in_triad=[ind_in_triad;iis_triad(line_number,:)];
            judge_i_triad=[];%存放新感染节点在iss_triad中所在行的三角形
            [line_number2,column_numbe2]=find(iss_triad==ind);
            judge_i_triad=[judge_i_triad;iss_triad(line_number2,:)];
            judge_triad=[];%存放新感染节点在triad中所在行的三角形
            [line_number3,column_numbe3]=find(sss_triad==ind);
            judge_triad=[judge_triad;sss_triad(line_number3,:)];
            if (susp(ind)==1)
                %判断新感染的节点(ind)是否是iis_triad中的易感节点，如果是，从iis_triad中删除这个三角形
                if (isempty(ind_in_triad)==0)
                    iii_triad=[iii_triad;ind_in_triad];
                    iis_triad(find(ismember(iis_triad,ind_in_triad,'rows')),:)=[];
                end
                %判断新感染的节点（ind）是否在iss_triad中
                if (isempty(judge_i_triad)==0)
                    iis_triad=[iis_triad;judge_i_triad];
                    iss_triad(find(ismember(iss_triad,judge_i_triad,'rows')),:)=[];
                end
                %判断新感染的节点是否在sss_triad中
                if (isempty(judge_triad)==0)
                    iss_triad=[iss_triad;judge_triad];
                    sss_triad(find(ismember(sss_triad,judge_triad,'rows')),:)=[];
                end
                susp(ind)=0;
                nodes = [nodes; ind];
            end

            %-----------------------------------------------------------------------------------------------------------------------------------------------------------
        elseif (next_inf2 < next_inf && next_inf2 < next_rec)
            time=next_inf2;
            IIS_I_ind=[];
            IIS_S_ind=[];
            for i=1:length(triad_inf)
                IIS_I_ind=[IIS_I_ind;nodes(find(nodes==triad_inf(i)))];%找到能够传染的三角形当中的两个传染者
            end
            IIS_S_ind=setdiff(triad_inf,IIS_I_ind);%找出三角形当中的S
            susp(IIS_S_ind)=0;
            nodes = [nodes; IIS_S_ind];
            iii_triad=[iii_triad;triad_inf];
            iis_triad(find(ismember(iis_triad,triad_inf,'rows')),:)=[];
            %判断易感节点是否在iis_triad中
            judge_iis_triad5=[];
            [line_number7,column_number7]=find(iis_triad==IIS_S_ind);
            judge_iis_triad5=[judge_iis_triad5;iis_triad(line_number7,:)];
            if (isempty(judge_iis_triad5)==0)
                iii_triad=[iii_triad;judge_iis_triad5];
                iis_triad(find(ismember(iis_triad,judge_iis_triad5,'rows')),:)=[];
            end
            %判断易感节点在不在i_triad中
            judge_i_triad5=[];
            [line_number8,column_number8]=find(iss_triad==IIS_S_ind);
            judge_i_triad5=[judge_i_triad5;iss_triad(line_number8,:)];
            if (isempty(judge_i_triad5)==0)
                iis_triad=[iis_triad;judge_i_triad5];
                iss_triad(find(ismember(iss_triad,judge_i_triad5,'rows')),:)=[];
            end
            %判断易感节点在不在triad中
            judge_triad5=[];
            [line_number9,column_number9]=find(sss_triad==IIS_S_ind);
            judge_triad5=[judge_triad5;sss_triad(line_number9,:)];
            if (isempty(judge_triad5)==0)
                iss_triad=[iss_triad;judge_triad5];
                sss_triad(find(ismember(sss_triad,judge_triad5,'rows')),:)=[];
            end


            %-----------------------------------------------------------------------------------------------------------------------------------------------------------
        elseif (next_rec < next_inf &&  next_rec < next_inf2)
            judge_iii_triad=[];
            [line_number4,column_number4]=find(iii_triad==node_rec);
            judge_iii_triad=[judge_iii_triad;iii_triad(line_number4,:)];
            judge_iis_triad=[];
            [line_number5,column_number5]=find(iis_triad==node_rec);
            judge_iis_triad=[judge_iis_triad;iis_triad(line_number5,:)];
            judge_i_triad2=[];
            [line_number6,column_number6]=find(iss_triad==node_rec);
            judge_i_triad2=[judge_i_triad2;iss_triad(line_number6,:)];
            if (susp(node_rec)==0)
                %判断恢复节点(node_rec)在不在iii_triad中
                if (isempty(judge_iii_triad)==0)
                    iii_triad(find(ismember(iii_triad,judge_iii_triad,'rows')),:)=[];
                    iis_triad=[iis_triad;judge_iii_triad];
                end
                %判断恢复节点在不在iis_triad中
                if (isempty(judge_iis_triad)==0)
                    iss_triad=[iss_triad;judge_iis_triad];
                    iis_triad(find(ismember(iis_triad,judge_iis_triad,'rows')),:)=[];
                end
                %判断恢复节点在不在i_triad中
                if (isempty(judge_i_triad2)==0)
                    sss_triad=[sss_triad;judge_i_triad2];
                    iss_triad(find(ismember(iss_triad,judge_i_triad2,'rows')),:)=[];
                end
                susp(node_rec)=1;
                nodes(ind_rec)=[];
                Tinf=[];
                Tinf2=[];
                Trec=[];
                time=next_rec;
            end
        end

        time

        %------------------------找每个时刻的[SI],[II],[I]的值--------------------------
        s_ind=find(susp==1);
        i_ind=find(susp==0);
        SS=full(sum(sum(A1(s_ind,s_ind))));
        SI=full(sum(sum(A1(s_ind,i_ind))));
        II=full(sum(sum(A1(i_ind,i_ind))));
        I=length(i_ind);
        IIS=size(iis_triad,1);
        ee=[ee;time,SS];%第1列是时间，第2列是每个时刻对应的SS的值
        cc=[cc; time,SI];
        gg=[gg; time,II];
        xx=[xx; time,I];
        mm=[mm;time,IIS];
 
    end
    Infect_I_50_1(m)=xx(end,2);
    fprintf('beta=%f时的I已保存！\n',beta);
end
plot(Beta1,Infect_I_50_1,'b*','LineWidth',3);
%plot(xx(:,1),xx(:,2),'color',[.5 .5 .5]);          %I       xx(:,1)是时间；xx(:,2)是I
save dataChaosSimSIS_I_50_1

% xlabel('t');
% ylabel('I');
% xlim([0,500]);
% ylim([0,2000]);
% set(gca,'yticklabel',{'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1'});
% set(gca,'color','none');%设置背景色为透明
% hold on
% h1=plot([0,500],[400,400],'--k'); %两点画直线
% h2=plot([0,500],[1000,1000],'--r');
% legend([h1,h2],'\rho_2_-^*','\rho_2_+^*','Location','NorthEast');
% % legend('N=1000,k_1=6,k_2=3,\beta=0.5,\beta_2=0.1,\gamma=0.05,I(0)=50','Location','SouthEast');
   
% hold on
% figure(2)
% plot(mm(:,1),mm(:,2),'color',[.5 .5 .5]);            %IIS
% xlabel('t');
% ylabel('IIS');
% legend('N=1000,k_1=6,k_2=3k_1=3,k_2=3,\beta=0.5,\beta_2=0.1,\gamma=0.05,I(0)=50','Location','NorthEast');
%
% hold on
% figure(3)
% plot(cc(:,1),cc(:,2),'color',[.5 .5 .5]);            %SI
% xlabel('t');
% ylabel('SI');
% legend('N=1000,k_1=6,k_2=3,\beta=0.5,\beta_2=0.1,\gamma=0.05,I(0)=50','Location','NorthEast');
%
% hold on
% figure(4)
% plot(gg(:,1),gg(:,2),'color',[.5 .5 .5]);   %II
% xlabel('t');
% ylabel('II');
% legend('N=1000,k_1=6,k_2=3,\beta=0.5,\beta_2=0.1,\gamma=0.05,I(0)=50','Location','SouthEast');
%
% hold on
% figure(5)
% plot(ee(:,1),ee(:,2),'color',[.5 .5 .5]);   %II
% xlabel('t');
% ylabel('SS');
% legend('N=1000,k_1=6,k_2=3,\beta=0.5,\beta_2=0.1,\gamma=0.05,I(0)=50','Location','NorthEast');
%
% save dataChaosSimSIS1
