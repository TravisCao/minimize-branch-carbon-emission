clc;
clear;

define_constants;

loadProfile = load("data/Case39_LoadData_120Hours.mat");
HourlyLoadProfile = mean(reshape(loadProfile.PD, [39, 24, 5]), 3);

HourlyLoadProfile = HourlyLoadProfile(:, 1:2:24);

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


T = size(HourlyLoadProfile, 2);
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
% 上网电价 常数 2
renew_gencost_general = [2 0 0 3 0 2 0];
mpc_modified.gen = [wind_gen; solar_gen; mpc_modified.gen];

% adjust gencost original, 10 is original gens number
mpc_modified.gencost = repmat([2	0	0	3	0.001	0.2	0.2], 10, 1);
mpc_modified.gencost = [renew_gencost_general; renew_gencost_general; mpc_modified.gencost];
mpc_modified.gencost(renew.hydro_gen, :) = repmat(renew_gencost_general, length(renew.hydro_gen), 1);

% 新能源
Eg_modified = Eg_orig;
Eg_modified = [0; 0; Eg_modified];
Eg_modified(renew.hydro_gen) = 0;

emission_opf = CEF_case39(T, mpc_modified, Eg_modified, HourlyLoadProfile, renew, 'opf', []);

sprintf("emission OPF is %.3f", mean(emission_opf))
% TODO 计算cost & 新能源出力情况

% "emission OPF is 380.609"

% ------- 新能源完全不调 -- 设定新能源max -- 负负荷 用OPF

% not negative load 

% % 去除机组设定
% mpc_renew_max = mpc;
% mpc_renew_max.gen = mpc_renew_max.gen(4:end, :);
% mpc_renew_max.gencost = mpc_renew_max.gencost(4:end, :);
% 
% 
% % 似乎不需要考虑去掉 slack node 的问题
% % mpc_renew_max.bus(31, 2) = 2;
% % mpc_renew_max.bus(30, 2) = 3;
% 
% Eg_renew_max = Eg_orig;
% Eg_renew_max = Eg_renew_max([1, 4:9]);

nGens = size(mpc_modified.gen, 1);

mpc_renew_max = mpc_modified;

% renewable no cost
mpc_renew_max.gencost([1 2 4 5 12], :) = repmat([2 0 0 3 0 0 0], 5, 1);

% Pg_renew_max = ones(nGens * T, 1);
% for t = 1:T
%     sum_renewable = solarProfile(t) + hydroProfile(t) * 3 + windProfile(t);
%     Pg_renew_max((t - 1) * nGens + 1: t * nGens) = (sum(HourlyLoadProfile(:, t)) - sum_renewable) / 7;
%     Pg_renew_max((t - 1) * nGens + renew.solar_gen) = solarProfile(t);
%     Pg_renew_max((t - 1) * nGens + renew.hydro_gen) = hydroProfile(t);
%     Pg_renew_max((t - 1) * nGens + renew.wind_gen) = windProfile(t);
% end

emission_renew_max = CEF_case39(T, mpc_renew_max, Eg_modified, HourlyLoadProfile, renew, 'max-renew', []);

sprintf("emission renew max is %.3f", mean(emission_renew_max))

% TODO 计算cost & 新能源出力情况

% "emission renew max is 129.365"

% TODO

% ------- 低碳调度

% mpc_modified.bus(31, 2) = 2;
% mpc_modified.bus(1, 2) = 3;

% 使用和opf一样的 mpc_modified（修改3个hydro gen & 添加 1 wind和 1 solar）
nGens = size(mpc_modified.gen, 1);
nBranches = size(mpc_modified.branch, 1);
nvars = T * nGens; % T * 机组数量
ObjectiveFunction = @(Pg) CEF_case39(T, mpc_modified, Eg_modified, HourlyLoadProfile, renew, 'pf', Pg);

PTDF = makePTDF(mpc_modified);

A = zeros(T * nBranches, nvars);
b = zeros(T * nBranches, 1);

for t = 1:T
    A((t - 1) * nBranches + 1 : t * nBranches, (t - 1) * nGens + 1 : t * nGens) = PTDF(:, mpc_modified.gen(:, GEN_BUS));
    b((t - 1) * nBranches + 1 : t * nBranches) = PTDF * HourlyLoadProfile(:, t) + mpc.branch(:, RATE_A);
end

% 初始化 Aeq 和 beq
Aeq = zeros(T, nvars);
beq = zeros(T, 1);

% 构造 Aeq 和 beq
for t = 1:T
    Aeq(t, (t - 1) * nGens + 1 : t * nGens) = ones(1, nGens);
    beq(t) = sum(HourlyLoadProfile(:, t));
end

% lower bound & upper bound
lb = zeros(nvars, 1);
ub = repmat(mpc_modified.gen(:, PMAX), T, 1);
for t = 1:T
%     ub((t - 1) * nGens + 1: t * nGens) = 900;
    ub((t - 1) * nGens + renew.solar_gen) = solarProfile(t);
    ub((t - 1) * nGens + renew.hydro_gen) = hydroProfile(t);
    ub((t - 1) * nGens + renew.wind_gen) = windProfile(t);
end

meanPg = sum(mpc.bus(:, 3)) / size(mpc.gen, 1);
rng(42) % For reproducibility
options = optimoptions("ga",'ConstraintTolerance', 1e-5, ...
    'PlotFcn', {'gaplotbestf','gaplotrange'}, 'Display', 'diagnose', ...
    'MaxGenerations', 800);

% 初始设置每个新能源机组为0，其他机组发电机组为800

initialPopulationMatrix = ones(nvars, 1)'* 550;
for t = 1:T
    sum_renewable = (solarProfile(t) + hydroProfile(t) * 3 + windProfile(t));
    initialPopulationMatrix((t - 1) * nGens + 1: t * nGens) = (sum(HourlyLoadProfile(:, t)) - sum_renewable * 1/2) / 7;
    initialPopulationMatrix((t - 1) * nGens + renew.solar_gen) = solarProfile(t) * 1/2;
    initialPopulationMatrix((t - 1) * nGens + renew.hydro_gen) = hydroProfile(t) * 1/2;
    initialPopulationMatrix((t - 1) * nGens + renew.wind_gen) = windProfile(t) * 1/2;
end

% initialPopulationRange = zeros(nvars, 2)';
% for t = 1:T
%     initialPopulationRange(1, (t - 1) * nGens + 1: t * nGens) = 500;
%     initialPopulationRange(2, (t - 1) * nGens + 1: t * nGens) = 900;
%     initialPopulationRange(1, (t - 1) * nGens + renew.solar_gen) = solarProfile(t) * 1/2;
%     initialPopulationRange(2, (t - 1) * nGens + renew.solar_gen) = solarProfile(t);
%     initialPopulationRange(1, (t - 1) * nGens + renew.hydro_gen) = hydroProfile(t) * 1/2;
%     initialPopulationRange(2, (t - 1) * nGens + renew.hydro_gen) = hydroProfile(t) * 1;
%     initialPopulationRange(1, (t - 1) * nGens + renew.wind_gen) = windProfile(t) * 1/2;
%     initialPopulationRange(2, (t - 1) * nGens + renew.wind_gen) = windProfile(t);
% end

options.InitialPopulationMatrix = initialPopulationMatrix;
options.InitialPopulationRange = initialPopulationRange;
options.CrossoverFraction = 0.8;
options.PopulationSize = 400;
[x, fval, exitflag, output, population, scores] = ga(ObjectiveFunction, nvars, A, b, Aeq, beq, lb, ub, [], options)