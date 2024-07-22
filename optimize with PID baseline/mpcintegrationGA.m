% 清除工作区变量
clear; clc;

% 设置固定的随机种子
rng(62); % 42是种子值，可以选择任何整数

%% 定义系统参数
mp = 0.024;
mr = 0.095;
Lp = 0.129;
Lr = 0.085;
Jp = 3.33 * 10^(-5);
Jr = 5.72 * 10^(-5);
kt = 0.042;
km = 0.042;
Rm = 8.4;
g = 9.8;
JT = Jp * mp * Lr^2 + Jr * Jp + 1/4 * (Jr * mp * Lp^2);

A_c = [0 0 1 0;
       0 0 0 1;
       0 mp^2 * Lp^2 * Lr * g / (4 * JT) 0 0;
       0 mp * g * Lp * (Jr + mp * Lr^2) / (2 * JT) 0 0];

B_c = [0; 0; (Jp + 1/4 * (mp * Lp^2)) / JT; mp * Lp * Lr / (2 * JT)];
B_c = (kt / Rm) * B_c;
A_c(3, 3) = A_c(3, 3) - (kt * km / Rm);
A_c(4, 3) = A_c(4, 3) - (kt * km / Rm);
C = eye(4); % 观测所有状态
D = zeros(4, 1);

%% 离散化系统
Ts = 0.05; % 调整采样时间
sys_c = ss(A_c, B_c, C, D);
sys_d = c2d(sys_c, Ts);

A = sys_d.A;
B = sys_d.B;

% 检查离散化后的系统是否稳定
eig_values = eig(A);
if all(abs(eig_values) < 1)
    disp('系统稳定');
else
    disp('系统不稳定');
end


%% GA优化设置
pid_baseline = 0.0473; % 控制输入积分目标
options = optimoptions('ga', 'PopulationSize', 20, 'MaxGenerations', 50, 'Display', 'iter', 'UseParallel', false); % 禁用并行计算
lb = [0, 0, 0, 0, 0, 0]; % 下边界
ub = [100, 100, 100, 100, 1, 1]; % 上边界

% 目标函数定义
objectiveFunction = @(params) objectiveFunctionGA(params, A, B, Ts, pid_baseline);

% 设置随机种子
rng(52);

% GA优化
[params_opt, fval] = ga(objectiveFunction, 6, [], [], [], [], lb, ub, [], options);

%% 保存优化结果

% 提取优化后的权重
Q_opt = diag(params_opt(1:4));
R_opt = params_opt(5);
R_delta_opt = params_opt(6);

% 打印优化后的权重矩阵
fprintf('Optimized Q Matrix:\n');
disp(Q_opt);
fprintf('Optimized R:\n');
disp(R_opt);
fprintf('Optimized R_delta:\n');
disp(R_delta_opt);

%% 创建MPC控制器对象
prediction_horizon = 20; % 调整预测时域 (Np)
control_horizon = 5; % 调整控制时域 (Nu)

mpc_obj = mpc(sys_d, Ts, prediction_horizon, control_horizon);

% 设置控制输入约束
mpc_obj.MV.Min = -10; % 控制输入（电机电压）的最小值
mpc_obj.MV.Max = 10;

% 应用优化后的权重矩阵
mpc_obj.Weights.OutputVariables = params_opt(1:4);
mpc_obj.Weights.ManipulatedVariables = R_opt;
mpc_obj.Weights.ManipulatedVariablesRate = R_delta_opt;

%% 初始状态和仿真设置
x0 = [1 * pi / 180; 1 * pi / 180; 0; 0]; % 初始状态
k_steps = 200; % 增加仿真步数
T_sim = k_steps * Ts; % 仿真时间

% 创建MPC状态对象
mpcstate_obj = mpcstate(mpc_obj);

% 参考信号
r = zeros(4, 1);

% 仿真数据
X_K = zeros(size(A, 1), k_steps + 1);
U_K = zeros(size(B, 2), k_steps);
X_K(:, 1) = x0;

% 仿真MPC控制器
for k = 1:k_steps
    % 计算MPC控制输入
    u = mpcmove(mpc_obj, mpcstate_obj, X_K(:, k), r);
    U_K(:, k) = u;
    
    % 更新状态
    X_K(:, k + 1) = A * X_K(:, k) + B * u;
end

% 计算MPC控制输入的总和
mpc_control_input_sum = sum(abs(U_K)) * Ts;

% 打印MPC控制输入的总和和PID基线
fprintf('MPC Control Input Sum: %.4f\n', mpc_control_input_sum);
fprintf('PID Control Input Baseline: %.4f\n', pid_baseline);

% 绘制系统响应图
time = (0:k_steps) * Ts; % 时间向量
figure;

subplot(2, 1, 1);
hold on;
for i = 1 : size(X_K, 1)
    plot(time, X_K(i, :));
end
title('Closed-Loop System States Response with MPC');
xlabel('Time (s)');
ylabel('States');
legend('\theta (Arm Angle)', '\alpha (Pendulum Angle)', ' (Arm Angular Velocity)', '(Pendulum Angular Velocity)');
hold off;

subplot(2, 1, 2);
hold on;
for i = 1 : size(U_K, 1)
    plot(time(1:end-1), U_K(i, :));
end
title('Control Input');
xlabel('Time (s)');
ylabel('u');
legend('u');
hold off;

% 绘制零极点图
figure;
pzmap(sys_d);
title('Discrete-Time System Pole-Zero Map');

% 输出特征值
disp('特征值:');
disp(eig_values);

% 计算各个状态的最大值和最小值
state_max = max(X_K, [], 2);
state_min = min(X_K, [], 2);

% 打印结果
fprintf('Maximum and Minimum values of each state:\n');
for i = 1:size(X_K, 1)
    fprintf('State %d: Max = %.4f, Min = %.4f\n', i, state_max(i), state_min(i));
end

% 计算控制输入的最大值和最小值
u_max = max(U_K, [], 2);
u_min = min(U_K, [], 2);

% 打印结果
fprintf('Maximum and Minimum values of each control input:\n');
for i = 1:size(U_K, 1)
    fprintf('U %d: Max = %.4f, Min = %.4f\n', i, u_max(i), u_min(i));
end

% 目标函数定义
function J = objectiveFunctionGA(params, A, B, Ts, target_integral_abs)
    % 提取权重矩阵参数
    Q = diag(params(1:4));
    R = params(5);
    R_delta = params(6);

    % 创建MPC控制器对象
    prediction_horizon = 15; % 预测时域 (Np)
    control_horizon = 10; % 控制时域 (Nu)
    mpc_obj = mpc(ss(A, B, eye(size(A, 1)), zeros(size(A, 1), size(B, 2))), Ts, prediction_horizon, control_horizon);

    % 设置控制输入约束
    mpc_obj.MV.Min = -10; % 控制输入（电机电压）的最小值
    mpc_obj.MV.Max = 10;

    % 应用权重矩阵
    mpc_obj.Weights.OutputVariables = params(1:4);
    mpc_obj.Weights.ManipulatedVariables = R;
    mpc_obj.Weights.ManipulatedVariablesRate = R_delta;

    % 初始状态
    x0 = [1 * pi / 180; 1 * pi / 180; 0; 0]; % 初始状态
    k_steps = 100; % 初步优化时减少仿真步数

    % 创建MPC状态对象
    mpcstate_obj = mpcstate(mpc_obj);

    % 参考信号
    r = zeros(size(A, 1), 1);

    % 仿真数据
    X_K = zeros(size(A, 1), k_steps + 1);
    U_K = zeros(size(B, 2), k_steps);
    X_K(:, 1) = x0;

    % 仿真MPC控制器
    for k = 1:k_steps
        % 计算MPC控制输入
        u = mpcmove(mpc_obj, mpcstate_obj, X_K(:, k), r);
        U_K(:, k) = u;
        
        % 更新状态
        X_K(:, k + 1) = A * X_K(:, k) + B * u;
    end

    % 计算控制输入的总和
    control_input_sum = sum(abs(U_K)) * Ts;

    % 计算目标函数值（绝对误差）
    J = abs(control_input_sum - target_integral_abs);
end
