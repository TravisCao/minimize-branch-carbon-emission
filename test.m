clc;
clear;

define_constants;

% a = 100;
% b = 1; % define constant values
% lb = [-3,-3];
% ub = [3,3];
% numberOfVariables = 2;
% FitnessFunction = @(x) parameterized_fitness(x,a,b);
% [x,fval] = ga(FitnessFunction,numberOfVariables,[],[],[],[],lb,ub);


T = 1;    % 持续时间 h
mpc = loadcase('case39');
n_bus = size(mpc.bus, 1);
Gn = length(mpc.gen(:, 1));  %bsn = length(mpc.bus(:,1));   brn = length(mpc.branch(:,1));
Eg = zeros(Gn, 1);           %  ECEF intensity   %发电机碳排放强度
Eg(mpc.gen(:, 2) > 600) = 0.875;      % 单位：tCO2/M·Wh    
Eg(300 < mpc.gen(:,2) & mpc.gen(:, 2) <= 600) = 0.500;  

nvars = T * size(mpc.gen, 1); % T * 机组数量
ObjectiveFunction = @(Pg) CEF_case39(T, mpc, Eg, Pg);
% ConstraintFunction = @(Pg) para_costr(mpc, Pg, mpc.bus(:, PD));

PTDF = makePTDF(mpc);

A = PTDF(:, end - 9:end);
b = PTDF * mpc.bus(:, PD) + mpc.branch(:, RATE_A);
Aeq = ones(nvars, 1)';
beq = sum(mpc.bus(:, PD));

% lower bound & upper bound
lb = mpc.gen(:, PMIN);
ub = mpc.gen(:, PMAX);

meanPg = sum(mpc.bus(:, 3)) / size(mpc.gen, 1);
rng default % For reproducibility
options = optimoptions("ga",'PlotFcn', {@gaplotbestf,@gaplotmaxconstr}, 'Display', 'iter');
options.InitialPopulationMatrix = ones(size(mpc.gen(:, 2)))' * sum(mpc.bus(:, 3)) / size(mpc.gen, 1);
options.PopulationSize = 50;
[x, fval] = ga(ObjectiveFunction, nvars, A, b, Aeq, beq, lb, ub, [], options)

