% 清除工作区变量
clear; clc;
%试过了与种子无关
%rng(62);

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
Vm = 18;

% 定义状态空间模型
A = [0 0 1 0;
     0 0 0 1;
     0 mp^2*Lp^2*Lr*g/(4*JT) 0 0;
     0 mp*g*Lp*(Jr+mp*Lr^2)/(2*JT) 0 0];
B = [0; 0; (Jp + 1/4*(mp*Lp^2))/JT; mp*Lp*Lr/(2*JT)];
C = eye(4);
D = zeros(4, 1);

% 检查系统的可控性
C_ctrb = [B A*B A*A*B A*A*A*B];
rank_CM = rank(C_ctrb);

if (rank_CM == size(A, 1))
    disp('The controllability matrix has full rank and the system is therefore controllable.');
else
    disp('The controllability matrix does not have full rank and the system is therefore not controllable.');
end

% 改变参数，将控制输入 T 转换为 V_m
B = (kt/Rm) * B;
A(3, 3) = A(3, 3) - (kt*km/Rm);
A(4, 3) = A(4, 3) - (kt*km/Rm);

% PSO优化设置
Ts = 0.01; % 时间步长
T = 0:Ts:10; % 仿真时间
x0 = [1 * pi / 180; 1 * pi / 180; 0; 0]; % 初始条件
target_integral_abs = 0.0473; % 目标控制输入积分绝对值

% 定义目标函数
objectiveFunction = @(params) objectiveFunctionPSO(params, A, B, Ts, target_integral_abs, T, x0);

% 初始猜测值（极点位置），包括复数部分
initialParams = [-30 + 10i, -30 - 10i, -10 + 5i, -10 - 5i];

% PSO优化
options = optimoptions('particleswarm', 'SwarmSize', 30, 'MaxIterations', 300, 'Display', 'iter');
lb = [-100, -100, -100, -100, -50, -50, -50, -50]; % 下边界，实部和虚部分开处理
ub = [0, 0, 0, 0, 50, 50, 50, 50]; % 上边界，实部和虚部分开处理

[params_opt, fval] = particleswarm(objectiveFunction, 8, lb, ub, options);

% 使用优化后的极点计算状态反馈增益
params_opt_complex = params_opt(1:4) + 1i * params_opt(5:8);
K_new = place(A, B, params_opt_complex);
K_real = real (K_new);
% 显示优化后的极点位置
disp('Optimized Pole Locations:');
disp(params_opt_complex);

% 定义闭环系统
A_cl = A - B * K_new;
sys_cl = ss(A_cl, B, C, D);

% 仿真优化后的极点配置控制器
[y, t, x] = initial(sys_cl, x0, T);
u = -K_new * x';

% 计算优化后的控制输入积分
integral_opt = sum(abs(u)) * Ts;
fprintf('Optimized Control Input Integral: %.4f\n', integral_opt);

% 计算各个状态的最大值和最小值
state_max = max(x);
state_min = min(x);

% 打印结果
fprintf('Maximum and Minimum values of each state:\n');
for i = 1:size(x, 2)
    fprintf('State %d: Max = %.4f, Min = %.4f\n', i, state_max(i), state_min(i));
end

% 绘制响应图
figure;
subplot(2, 1, 1);
plot(t, x);
title('Closed-Loop System States Response with Optimized Pole Placement');
xlabel('Time (s)');
ylabel('States');
legend('\theta (Arm Angle)', '\alpha (Pendulum Angle)', ' (Arm Angular Velocity)', '(Pendulum Angular Velocity)');

subplot(2, 1, 2);
plot(t, u);
title('Control Input');
xlabel('Time (s)');
ylabel('u');

% 绘制零极点图
figure;
pzmap(sys_cl);
title('Closed-Loop System Pole-Zero Map with Optimized Pole Placement');

% 目标函数定义
function J = objectiveFunctionPSO(params, A, B, Ts, target_integral_abs, T, x0)
    % 使用传递函数得到状态反馈增益
    try
        params_complex = params(1:4) + 1i * params(5:8);
        K = place(A, B, params_complex);
    catch
        J = Inf;
        return;
    end
    
    % 闭环系统矩阵
    A_cl = A - B * K;

    % 定义闭环系统
    sys_cl = ss(A_cl, B, eye(4), zeros(4, 1));

    % 仿真闭环系统
    [~, ~, x] = initial(sys_cl, x0, T);
    u = -K * x';

    % 计算控制输入的积分
    integral_u = sum(abs(u)) * Ts;

    % 目标函数：控制输入积分与目标值的绝对差
    J = abs(integral_u - target_integral_abs);
end
