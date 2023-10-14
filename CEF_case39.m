% function [En, Rg, RL, Rb1, Rl, CE_Gen, CE_Load, CE_Loss, obj] = CEF_case39(Pg)
function [obj] = CEF_case39(T, mpc, Eg, HourlyLoadProfile, renew, choice, Pg)
    %……输入变量 
    %T—时间间隔，mpc—电网参数，Eg—发电机排放强度
    %……输出变量
    %En—ICEF intensity节点碳排放强度,Rg—ECEF机组碳排放rate
    %RL—load ICEF rate,Rb1—BCEF线路rate,Rl—BCEL线损rate
    %CE_Gen—机组排放量,CE_Load—负荷排放量,CE_Loss-线损排放总量  CE_Gen-CE_Load=CE_Loss；

    % 目标最小化区间 [1, 39], [4, 5], [13 14]
    targetBranch_row = [1 4 13];
    targetBranch_col = [39 5 14];
    sz = [39, 39];
    target_ind = sub2ind(sz, targetBranch_row, targetBranch_col);

    windProfile = renew.wind;
    solarProfile = renew.solar;
    hydroProfile = renew.hydro;

    %% 主程序
    opt = mpoption('out.all', 0, 'verbose', 0);
    % results =  runopf(mpc, opt);
    obj = zeros(T, 1);
    nGens = size(mpc.gen, 1);
    unconvergenceCount = 0;
    for t = 1:T
        
        % logging
        sprintf("Time is %d", t)

        mpc.bus(:, 3) = HourlyLoadProfile(:, t);

        if strcmp(choice, 'max-renew')
            % set negative load as renewable max Pg
            mpc.gen(renew.wind_gen, 9) = windProfile(t);
            mpc.gen(renew.solar_gen, 9) = solarProfile(t);
            mpc.gen(renew.hydro_gen, 9) = hydroProfile(t);
%             mpc.bus(renew.hydro_bus, 3) = mpc.bus(renew.hydro_bus, 3) - hydroProfile(t);
%             mpc.bus(renew.wind_bus, 3) = mpc.bus(renew.wind_bus, 3) - windProfile(t);
%             mpc.bus(renew.solar_bus, 3) = mpc.bus(renew.solar_bus, 3) - solarProfile(t);
            [results, success] = rundcopf(mpc, opt);
            obj(t) = CEF_hourly(results, success, 'no');
            sprintf("run max renew")
            results.gen(:, 2)
        else
            % set renewable PMAX to time-varient
            mpc.gen(renew.wind_gen, 9) = windProfile(t);
            mpc.gen(renew.solar_gen, 9) = solarProfile(t);
            mpc.gen(renew.hydro_gen, 9) = hydroProfile(t);
            if strcmp(choice, 'pf')
                mpc.gen(:, 2) = Pg((t - 1) * nGens + 1: t * nGens);
                [results, success] = rundcpf(mpc, opt);
                obj(t) = CEF_hourly(results, success, 'fitness');
                results.gen(:, 2)
            elseif strcmp(choice, 'opf')
                [results, success] = rundcopf(mpc);
                obj(t) = CEF_hourly(results, success, 'no');
                sprintf("run opf")
                results.gen(:, 2)
            end
        end
        if success == 0
            unconvergenceCount = unconvergenceCount + 1;
        end
    end
    
    unconvergenceCount
    obj
    
    if strcmp(choice, 'pf')
        obj = mean(obj);
    end

    function obj_hourly = CEF_hourly(results, success, flag)
        % result
        % [results, success] = runpf(mpc);
        
%         pg_result = results.gen(:, 2);
%         pg_ub = results.gen(:, 9);
%         pg_lb = results.gen(:, 10);
        if success == 0 
            obj_hourly = 100000;
            return
%         elseif  ~all(pg_result < pg_ub | abs(pg_result - pg_ub) < 1e-3)
%             obj_hourly = sum(abs(pg_result - pg_ub));
%             return
%         elseif ~all(pg_result > pg_lb | abs(pg_result - pg_ub) < 1e-3)
%             obj_hourly = sum(abs(pg_result - pg_lb));
%             return
        end
    
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
        obj_hourly = sum(Rb1(target_ind));
%         renew_gen = [1 2 4 5 12];
%         if strcmp(flag, 'fitness')
%             if abs(sum(mpc.gen(renew_gen, 9) - mpc.gen(renew_gen, 2))) > 100
%                 obj_hourly = obj_hourly + 1000 * sum(mpc.gen(renew_gen, 9) - mpc.gen(renew_gen, 2));
%             end
%         end
    end

end