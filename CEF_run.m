


%% load

% load in each timestep is calculated by using the average of each day in five days (120h)
% P_{i, t} = 1 / 5 * \sum_{d=1}^{5} P_{i, t + (d-1) * 24}

loadProfile = load("data/Case39_LoadData_120Hours.mat");
HourlyLoadProfile = mean(reshape(loadProfile.PD, [39, 24, 5]), 3);

windProfile = load("data/wind.mat").wind;
solarProfile = load("data/solar.mat").solar;
hydroProfile = load("data/hydro.mat").hydro;

% 这里设定的前提是将renew 加进gen，1:wind， 2:solar，12，4，5：hydro

renew.wind = windProfile;
renew.wind_bus = 3;
renew.wind_gen = 1;
renew.solar = solarProfile;
renew.solar_bus = 15;
renew.solar_gen = 2;
renew.hydro = hydroProfile;
renew.hydro_gen = [12 4 5];
renew.hydro_bus = [39 31 32];


%% renewable energy
% TODO

T = 24;
mpc = loadcase('case39');

% 增加其他机组的总发电量
% mpc.gen(5:end, 9) = 1800;
% 增加所有线路的总容量
% mpc.branch(:, 6) = 3000;

Gn = length(mpc.gen(:, 1));  %bsn = length(mpc.bus(:, 1));   brn = length(mpc.branch(:, 1));
Eg_orig = zeros(Gn, 1);           %  ECEF intensity   %发电机碳排放强度
Eg_orig(mpc.gen(:, 2) > 600) = 0.875;      % 单位：tCO2/M·Wh    
Eg_orig(300 < mpc.gen(:,2) & mpc.gen(:, 2) <= 600) = 0.500;

% 增加其他机组容量到1000
mpc.gen(:, PMAX) = 1100;
mpc.branch(:, RATE_A) = 4000;

% ------- 常规经济调度

mpc_modified = mpc;

% 加新能源
wind_gen = mpc.gen(1, :);
wind_gen(GEN_BUS) = 3;

solar_gen = mpc.gen(1, :);
solar_gen(GEN_BUS) = 15;

% 1 & 2 新加：风光机组, 12 4 5 丽水内部 水电
% 上网电价 常数 0.8
renew_gencost_general = [2 0 0 3 0 0 0.8];
mpc_modified.gen = [wind_gen; solar_gen; mpc_modified.gen];
mpc_modified.gencost = [renew_gencost_general; renew_gencost_general; mpc_modified.gencost];
mpc_modified.gencost(renew.hydro_gen, :) = repmat(renew_gencost_general, length(renew.hydro_gen), 1);

% 新能源
Eg_modified = Eg_orig;
Eg_modified = [0; 0; Eg_modified];
Eg_modified(renew.hydro_gen) = 0;

nGens = size(mpc_modified.gen, 1);

initialPopulationMatrix = ones(nvars, 1)';
for t = 1:T
    sum_renewable = solarProfile(t) + hydroProfile(t) * 3 + windProfile(t);
    initialPopulationMatrix((t - 1) * nGens + 1: t * nGens) = (sum(HourlyLoadProfile(:, t)) - sum_renewable) / 7;
    initialPopulationMatrix((t - 1) * nGens + renew.solar_gen) = solarProfile(t);
    initialPopulationMatrix((t - 1) * nGens + renew.hydro_gen) = hydroProfile(t);
    initialPopulationMatrix((t - 1) * nGens + renew.wind_gen) = windProfile(t);
end


mpc_modified.bus(:, 3) = HourlyLoadProfile(:, 1);
mpc_modified.gen(:, 2) = initialPopulationMatrix(1:12);

mpc_modified.gen(renew.wind_gen, 9) = windProfile(1);
mpc_modified.gen(renew.solar_gen, 9) = solarProfile(1);
mpc_modified.gen(renew.hydro_gen, 9) = hydroProfile(1);

[obj] = CEF(T, mpc_modified, Eg_modified)

mpc_renew_max.bus(:, 3) = HourlyLoadProfile(:, 1);
mpc_renew_max.bus(renew.hydro_bus, 3) = mpc_renew_max.bus(renew.hydro_bus, 3) - hydroProfile(1);
mpc_renew_max.bus(renew.wind_bus, 3) = mpc_renew_max.bus(renew.wind_bus, 3) - windProfile(1);
mpc_renew_max.bus(renew.solar_bus, 3) = mpc_renew_max.bus(renew.solar_bus, 3) - solarProfile(1);
sum_renewable = solarProfile(1) + hydroProfile(1) * 3 + windProfile(1);
mpc_renew_max.gen(:, 2) = (sum(HourlyLoadProfile(:, 1)) - sum_renewable) / 7;
CEF(T, mpc_renew_max, Eg_renew_max)



% % carbon emission rate 
% CR_g = Eg.*results.gen(:,2);  % tCO2/h
% CR_n = En.*results.bus(:,3);
% CR_loss = sum(sum(Rl));       % =sum(CR_g)-sum(CR_n)
% %carbon emission rate quantity
% CE_Gen = CR_g*T;
% CE_Load = CR_n*T;
% CE_Loss = CR_loss*T;
% sum(CE_Gen)
% %% 图形绘制
% %节点CEF intensity 
% figure
% b = bar(En); 
% b.FaceColor=[1 0.5 0];