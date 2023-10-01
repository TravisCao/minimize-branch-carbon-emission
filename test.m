clc;
clear;

define_constants;

% a = 100;
% b = 1; % define constant values
% lb = [-3,-3];
% ub = [3, 3];
% numberOfVariables = 2;
% FitnessFunction = @(x) parameterized_fitness(x,a,b);
% [x,fval] = ga(FitnessFunction, numberOfVariables,[],[],[],[],lb,ub);

mpc = loadcase('case39');
n_bus = size(mpc.bus, 1);
Gn = length(mpc.gen(:, 1));  %bsn = length(mpc.bus(:, 1));   brn = length(mpc.branch(:, 1));
Eg = zeros(Gn, 1);           %  ECEF intensity   %发电机碳排放强度
Eg(mpc.gen(:, 2) > 600) = 0.875;      % 单位：tCO2/M·Wh    
Eg(300 < mpc.gen(:,2) & mpc.gen(:, 2) <= 600) = 0.500;  

loadProfile = load("data/Case39_LoadData_120Hours.mat");
HourlyLoadProfile = mean(reshape(loadProfile.PD, [39, 24, 5]), 3);

T = 24;
nGens = size(mpc.gen, 1);
nvars = T * nGens; % T * 机组数量
ObjectiveFunction = @(Pg) CEF_case39(T, mpc, Eg, HourlyLoadProfile, Pg);
% ConstraintFunction = @(Pg) para_costr(mpc, Pg, mpc.bus(:, PD));

PTDF = makePTDF(mpc);

% 初始化A和b
A = zeros(T * size(PTDF, 1), nvars);
b = zeros(T * size(PTDF, 1), 1);

% 构造A和b
for t = 1:T
    A((t-1)*size(PTDF, 1) + 1 : t*size(PTDF, 1), (t-1)*nGens + 1 : t*nGens) = PTDF(:, end - nGens + 1:end);
    b((t-1)*size(PTDF, 1) + 1 : t*size(PTDF, 1)) = PTDF * HourlyLoadProfile(:, t) + mpc.branch(:, RATE_A);
end

% 初始化 Aeq 和 beq
Aeq = zeros(T, nvars);
beq = zeros(T, 1);

% 构造 Aeq 和 beq
for t = 1:T
    Aeq(t, (t-1)*nGens + 1 : t*nGens) = ones(1, nGens);
    beq(t) = sum(HourlyLoadProfile(:, t));
end

% lower bound & upper bound
% lb = repmat(mpc.gen(:, PMIN), T, 1);
% ub = repmat(mpc.gen(:, PMAX), T, 1);

lb = [];
ub = [];

meanPg = sum(mpc.bus(:, 3)) / size(mpc.gen, 1);
rng default % For reproducibility
options = optimoptions("ga",'PlotFcn', {@gaplotbestf,@gaplotmaxconstr}, 'Display', 'iter');

% 初始设置每个发电机组 500
options.InitialPopulationMatrix = repmat(ones(size(mpc.gen(:, 2))) * 500, T, 1)';
options.PopulationSize = 50;
[x, fval] = ga(ObjectiveFunction, nvars, A, b, Aeq, beq, lb, ub, [], options)

