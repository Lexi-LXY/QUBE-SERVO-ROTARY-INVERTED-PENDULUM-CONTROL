% 清除工作区变量
clear; clc;



% 定义系统参数
mp=0.024;
mr=0.095;
Lp=0.129;
Lr=0.085;
Jp=3.33*10^(-5);
Jr=5.72*10^(-5);
kt=0.042;
km=0.042;
Rm=8.4;
g=9.8;
JT=Jp*mp*Lr^2+Jr*Jp+1/4*(Jr*mp*Lp^2);

A_c = [0 0 1 0;
       0 0 0 1;
       0 mp^2*Lp^2*Lr*g/(4*JT) 0 0;
       0 mp*g*Lp*(Jr+mp*Lr^2)/(2*JT) 0 0];

B_c = [0; 0; (Jp+1/4*(mp*Lp^2))/JT; mp*Lp*Lr/(2*JT)];
B_c = (kt/Rm) * B_c;
A_c(3,3) = A_c(3,3) - (kt*km/Rm);
A_c(4,3) = A_c(4,3) - (kt*km/Rm);
C = eye(4); % only pendulum angle and arm angle
D = zeros(4, 1);

% 离散化系统
Ts = 0.05; % 采样时间
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

% 创建MPC控制器对象
prediction_horizon = 15; % %(Np) prediction horizon
control_horizon = 10; % %(Nu) control horizon

mpc_obj = mpc(sys_d, Ts, prediction_horizon, control_horizon);

mpc_obj.MV.Min = -10; % measured variable (i.e., motor voltage)
mpc_obj.MV.Max = 10; 
mpc_obj.Model.Plant.InputName = 'Motor_Voltage';
mpc_obj.Model.Plant.InputUnit = 'Voltage';
mpc_obj.Model.Nominal.X = [0;0;0;0];
mpc_obj.Model.Nominal.Y = [0;0;0;0];
mpc_obj.Model.Plant.OutputName = {'theta'; 'alpha';'thetadot';'alphadot'};
mpc_obj.Model.Plant.OutputUnit = {'degree'; 'degree';'degree/s';'degree/s'};
%这些设置定义了系统模型的输入名称、单位以及输出名称和单位。这些信息有助于MPC控制器正确地理解和处理系统的输入输出。

mpc_obj.OutputVariables(1).MinECR = 1;
mpc_obj.OutputVariables(1).MaxECR = 1;

mpc_obj.OutputVariables(2).Min = -0.05;
mpc_obj.OutputVariables(2).Max = 0.05;
mpc_obj.OutputVariables(2).MinECR = 1;
mpc_obj.OutputVariables(2).MaxECR = 1;

mpc_obj.Weights.ManipulatedVariables = 0.1;
mpc_obj.Weights.ManipulatedVariablesRate = 0.1;
mpc_obj.Weights.OutputVariables = [1 1 1 1];

%这些设置了输出变量的约束和权重：
%MinECR 和 MaxECR 定义了输出变量的容许误差限制。
%Weights.ManipulatedVariables 和 Weights.ManipulatedVariablesRate 设置了控制输入和其变化率的权重。
%Weights.OutputVariables 设置了输出变量的权重，这影响了控制器在优化过程中对输出变量的重视程度。

% 初始状态和仿真设置
x0 = [1 * pi / 180; 1 * pi / 180; 0; 0]; % 初始状态
k_steps = 200; % 仿真步数
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
    control_input_sum = sum(abs(U_K)) * Ts;

end
fprintf('MPC Control Input Sum: %.4f\n', control_input_sum);

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

% 打印权重矩阵
Q = mpc_obj.Weights.OutputVariables;
R = mpc_obj.Weights.ManipulatedVariables;
F = mpc_obj.Weights.ManipulatedVariablesRate;
 
fprintf('Q Matrix:\n');
disp(Q);

fprintf('R Matrix:\n');
disp(R);

fprintf('F Matrix:\n');
disp(F);

