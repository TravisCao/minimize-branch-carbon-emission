function [obj, En,Rg,RL,Rb1,Rl,CE_Gen,CE_Load,CE_Loss] = CEF(T,mpc,Eg)
    %……输入变量 
    %T—时间间隔，mpc—电网参数，Eg—发电机排放强度
    %……输出变量
    %En—ICEF intensity节点碳排放强度,Rg—ECEF机组碳排放rate
    %RL—load ICEF rate,Rb1—BCEF线路rate,Rl—BCEL线损rate
    %CE_Gen—机组排放量,CE_Load—负荷排放量,CE_Loss-线损排放总量  CE_Gen-CE_Load=CE_Loss；
    
    %% %%%%%%%%%%%%%
    % T = 1;    % 持续时间 h
    % mpc = loadcase('case30');
    
    % Gn = length(mpc.gen(:,1));  %bsn = length(mpc.bus(:,1));   brn = length(mpc.branch(:,1));
    % Eg = zeros(Gn,1);           %  ECEF intensity   %发电机碳排放强度
    % Eg(mpc.gen(:,2)>40) = 0.875;      % 单位：tCO2/M·Wh    
    % Eg(20<mpc.gen(:,2) & mpc.gen(:,2)<40) = 0.500;  


    % 目标最小化区间 [1, 39], [4, 5], [13 14]
    targetBranch_row = [1 4 13];
    targetBranch_col = [39 5 14];
    sz = [39, 39];
    target_ind = sub2ind(sz, targetBranch_row, targetBranch_col);
    
    %% 主程序
    results =  rundcpf(mpc);
    results.gen(:, 2)

    %% 节点分类  
    Ejn = results.gen(:,1);  gn = length(Ejn);
    bsn = length(results.bus(:,1));  brn = length(results.branch(:,1));
    Ijn = results.bus(:,1);
    Ijn(Ejn)=[];
    %% ECEF intensity 
    PG = zeros(bsn,gn); 
    for g=1:gn
        PG(Ejn(g),g) = results.gen(g,2);
    end
    Rg = PG*Eg;                % ECEF rate 
    %% ICEE intensity
    Br_o = results.branch;
    Pb1 = zeros(bsn);   %节点注入功率矩阵
    Pb = zeros(bsn);    %节点流出功率矩阵
    LoN = find(results.bus(:,3)>0);   ln = length(LoN);   
    Lcp = results.bus(LoN,3);
    PLo = zeros(bsn,ln); 
    for l=1:ln
        PLo(LoN(l),l) = Lcp(l);
    end
    Br_loss = zeros(brn,1);
    for l=1:brn
        if Br_o(l,14)>0 
            Pb1(Br_o(l,1),Br_o(l,2)) = abs(Br_o(l,16));
            Pb(Br_o(l,2),Br_o(l,1)) = Br_o(l,14);
        else
            Pb1(Br_o(l,2),Br_o(l,1)) = abs(Br_o(l,14));
            Pb(Br_o(l,1),Br_o(l,2)) = Br_o(l,16);
        end
    %     Br_loss(l) = abs(Br_o(l,14)+Br_o(l,16));
    end
    % aaaa = [Pb1;PG'];
    Pnd = diag(ones(1,bsn+gn)*[Pb1;PG']);
    Pnd_p = Pnd;   Pb1_p = Pb1;  Rg_p = Rg;  Pb_p = Pb;   PLo_p = PLo;
    if ~det(Pnd)
        dl = find(sum(Pnd)==0);
        Pnd_p(dl,:)=[];  Pnd_p(:,dl)=[];
        Pb1_p(dl,:)=[];  Pb1_p(:,dl)=[];
        Pb_p(dl,:)=[];   Pb_p(:,dl)=[];
        Rg_p(dl)=[];     PLo_p(dl,:) = [];
    end
    En_p = (inv(Pnd_p-Pb1_p'))*Rg_p;       % ICEE intensity
    %% ICEF rate  & BCEF rate  & BCEL rate
    RL_p = diag(En_p)*PLo_p;          % ICEF rate
    Rb1_p = diag(En_p)*Pb1_p;         % BCEF rate
    % aaaa = (Pb_p'-Pb1_p);
    Rl_p = diag(En_p)*(Pb_p'-Pb1_p);   % BCEL rate
    %% 还原删减矩阵
    dbl = zeros(1,ln);     dbu = zeros(1,bsn);
    RL=RL_p;  Rb1=Rb1_p;  Rl=Rl_p; En=En_p;
    
    if ~det(Pnd)
        for d=1:length(dl)
            RL = [RL(1:dl(d)-1,:);zeros(1,length(RL(1,:)));RL(dl(d):end,:)];
            %     RL_m = [RL_m(:,1:dl(d)-1),zeros(length(RL_m(:,1)),1),RL_m(:,dl(d):end)];
            Rb1 = [Rb1(1:dl(d)-1,:);zeros(1,length(Rb1(:,1)));Rb1(dl(d):end,:)];
            Rb1 = [Rb1(:,1:dl(d)-1),zeros(length(Rb1(:,1)),1),Rb1(:,dl(d):end)];
            Rl = [Rl(1:dl(d)-1,:);zeros(1,length(Rl(:,1)));Rl(dl(d):end,:)];
            Rl = [Rl(:,1:dl(d)-1),zeros(length(Rl(:,1)),1),Rl(:,dl(d):end)];
            En = [En(1:dl(d)-1);0;En(dl(d):end)];
        end
    end
    % carbon emission rate 
    CR_g = Eg.*results.gen(:,2);  % tCO2/h
    CR_n = En.*results.bus(:,3);
    CR_loss = sum(sum(Rl));       % =sum(CR_g)-sum(CR_n)
    % carbon emission rate quantity
    CE_Gen = CR_g*T;
    CE_Load = CR_n*T;
    CE_Loss = CR_loss*T;
    
    % target
    obj = sum(Rb1(target_ind));
end