% 清除工作区变量
clear; clc;
rng(62); % 设置固定的随机种子

% 定义系统参数
mp = 0.024;
mr = 0.095;
Lp = 0.129;
Lr = 0.085;
Jp = 3.33*10^(-5);
Jr = 5.72*10^(-5);
kt = 0.042;
km = 0.042;
Rm = 8.4;
g = 9.8;
JT = Jp*mp*Lr^2 + Jr*Jp + 1/4*(Jr*mp*Lp^2);

A = [0 0 1 0;
     0 0 0 1;
     0 mp^2*Lp^2*Lr*g/(4*JT) 0 0;
     0 mp*g*Lp*(Jr+mp*Lr^2)/(2*JT) 0 0];

B = [0; 0; (Jp + 1/4*(mp*Lp^2))/JT; mp*Lp*Lr/(2*JT)];
B = (kt/Rm) * B;
A(3,3) = A(3,3) - (kt*km/Rm);
A(4,3) = A(4,3) - (kt*km/Rm);
C=eye(4);
D=zeros(4,1);
n = size(A, 1);
p = size(B, 2);

%% 仿真设置
T_sim = 10; % 仿真总时间
dt = 0.01; % 仿真步长
t = 0:dt:T_sim; % 时间范围

X_K = zeros(n, length(t)); % 状态轨迹
initial_angle = 1 * pi / 180; % 初始状态 (1 度转换为弧度)
X_K(:, 1) = [initial_angle; initial_angle; 0; 0]; % 初始状态
U_K = zeros(p, length(t)); % 输入轨迹，输入数量为 p

%% PSO优化PID参数
nVars = 6; % PID参数的个数
lb = [-50, -10, -50, -50, -10, -50]; % 下界
ub = [50, 10, 50, 50, 10, 50]; % 上界

% PSO算法参数
options = optimoptions('particleswarm', ...
    'SwarmSize', 30, ... % 增加种群大小
    'MaxIterations', 300, ... % 增加最大迭代次数
    'Display', 'iter');

% 调用PSO进行优化
[optParams, fval] = particleswarm(@(params) pidObjective(params, A, B, C, D, dt, t), nVars, lb, ub, options);

% 提取优化后的PID参数
Kp_alpha = optParams(1);
Ki_alpha = optParams(2);
Kd_alpha = optParams(3);
Kp_theta = optParams(4);
Ki_theta = optParams(5);
Kd_theta = optParams(6);

fprintf('Optimized PID parameters:\n');
fprintf('Kp_alpha = %.2f, Ki_alpha = %.2f, Kd_alpha = %.2f\n', Kp_alpha, Ki_alpha, Kd_alpha);
fprintf('Kp_theta = %.2f, Ki_theta = %.2f, Kd_theta = %.2f\n', Kp_theta, Ki_theta, Kd_theta);

%% 仿真循环
x = [initial_angle; initial_angle; 0; 0]; % 初始状态
integral_error_alpha = 0;
previous_error_alpha = 0;
integral_error_theta = 0;
previous_error_theta = 0;
u_max = 10;
u_min = -10;
control_input_integral=0;
for k = 1:length(t)
    X_K(:, k) = x;
    error_alpha = 0 - x(1); % 摆角期望值为0
    error_theta = 0 - x(2); % 臂角期望值为0
    integral_error_alpha = integral_error_alpha + error_alpha * dt;
    integral_error_theta = integral_error_theta + error_theta * dt;
    derivative_error_alpha = (error_alpha - previous_error_alpha) / dt;
    derivative_error_theta = (error_theta - previous_error_theta) / dt;

    u_alpha = Kp_alpha * error_alpha + Ki_alpha * integral_error_alpha + Kd_alpha * derivative_error_alpha;
    u_theta = Kp_theta * error_theta + Ki_theta * integral_error_theta + Kd_theta * derivative_error_theta;
    u = u_alpha + u_theta;

    if u > u_max
        u = u_max;
    elseif u < u_min
        u = u_min;
    end

    U_K(k) = u;
    previous_error_alpha = error_alpha;
    previous_error_theta = error_theta;
    dx = A * x + B * u;
    x = x + dx * dt;
    % 计算控制输入的积分
    control_input_integral = control_input_integral + abs(u) * dt;
end

fprintf('PID Control Input Integral: %.4f\n', control_input_integral);






%%  说明本章中的脚本在simulink中执行的时候，对应的theta 和alpha参数是反过来的
%  说明本章中的脚本在simulink中执行的时候，对应的theta 和alpha参数是反过来的
%  说明本章中的脚本在simulink中执行的时候，对应的theta 和alpha参数是反过来的
%  重要的事情说三遍

%% 绘制系统响应图
figure;
subplot(2, 1, 1);
plot(t, X_K);
title('Closed-Loop System State Response with Optimized Multi-Variable PID Control');
xlabel('Time (s)');
ylabel('States');
legend('\theta (Arm Angle)', '\alpha (Pendulum Angle)', ' (Arm Angular Velocity)', '(Pendulum Angular Velocity)');

subplot(2, 1, 2);
plot(t, U_K); % 控制输入向量
title('Control Input');
xlabel('Time (s)');
ylabel('u');

%% 计算各个状态的最大值和最小值
state_max = max(X_K, [], 2);
state_min = min(X_K, [], 2);

% 打印结果
fprintf('Maximum and Minimum values of each state:\n');
for i = 1:n
    fprintf('State %d: Max = %.4f, Min = %.4f\n', i, state_max(i), state_min(i));
end


%%
%% PID优化目标函数
function J = pidObjective(params, A, B, C, D, dt, t)
    % 解包PID参数
    Kp_alpha = params(1);
    Ki_alpha = params(2);
    Kd_alpha = params(3);
    Kp_theta = params(4);
    Ki_theta = params(5);
    Kd_theta = params(6);

    % 仿真设置
    x = [1 * pi / 180; 1 * pi / 180; 0; 0]; % 初始状态
    integral_error_alpha = 0;
    previous_error_alpha = 0;
    integral_error_theta = 0;
    previous_error_theta = 0;
    U = zeros(1, length(t));
    X = zeros(length(x), length(t));

    % 控制输入限制
    u_max = 10;
    u_min = -10;

    for k = 1:length(t)
        X(:, k) = x;
        error_alpha = 0 - x(1); % 摆角期望值为0
        error_theta = 0 - x(2); % 臂角期望值为0
        integral_error_alpha = integral_error_alpha + error_alpha * dt;
        integral_error_theta = integral_error_theta + error_theta * dt;
        derivative_error_alpha = (error_alpha - previous_error_alpha) / dt;
        derivative_error_theta = (error_theta - previous_error_theta) / dt;

        % 计算控制输入 (电压，单位: 伏特)
        u_alpha = Kp_alpha * error_alpha + Ki_alpha * integral_error_alpha + Kd_alpha * derivative_error_alpha;
        u_theta = Kp_theta * error_theta + Ki_theta * integral_error_theta + Kd_theta * derivative_error_theta;
        u = u_alpha + u_theta;
        
        % 应用控制输入限制
        if u > u_max
            u = u_max;
        elseif u < u_min
            u = u_min;
        end
        
        U(k) = u; % 记录控制输入 (电压)
        previous_error_alpha = error_alpha;
        previous_error_theta = error_theta;
        dx = A * x + B * u;
        x = x + dx * dt;
    end

    % 计算性能指标（如IAE）
    J = sum(abs(X(1, :))) + sum(abs(X(2, :))) + 0.1 * sum(U.^2); % 综合摆角和臂角的误差及控制输入的能量
end



