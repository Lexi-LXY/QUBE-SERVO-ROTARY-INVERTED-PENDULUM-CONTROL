% 清除工作区变量
clear; clc;

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

if (rank_CM == size(C_ctrb, 1))
    disp('The controllability matrix has full rank and the system is therefore controllable.');
else
    disp('The controllability matrix does not have full rank and the system is therefore not controllable.');
end

% 改变参数，将控制输入 T 转换为 V_m
B = (kt/Rm) * B;
A(3, 3) = A(3, 3) - (kt*km/Rm);
A(4, 3) = A(4, 3) - (kt*km/Rm);

% 手动选择极点位置
P_new = [-30, -10 + 10i, -10 - 10i, -1];

% 使用 acker 函数计算状态反馈增益
K_new = acker(A, B, P_new);

% 验证极点是否满足特征多项式
char_poly = poly(P_new);
disp('Desired characteristic polynomial coefficients:');
disp(char_poly);

% 闭环系统矩阵
A_cl = A - B * K_new;
char_poly_cl = poly(eig(A_cl));
disp('Closed-loop characteristic polynomial coefficients:');
disp(char_poly_cl);

% 检查特征多项式是否匹配
if isequal(round(char_poly, 4), round(char_poly_cl, 4))
    disp('The chosen poles satisfy the characteristic polynomial.');
else
    disp('The chosen poles do NOT satisfy the characteristic polynomial.');
end

% 定义闭环系统
sys_cl = ss(A_cl, B, C, D);

% 绘制极点配置控制器的系统响应图
t = 0:0.01:10; % 时间范围
x0 = [1 * pi / 180; 1 * pi / 180; 0; 0]; % 初始条件
[y, t, x] = initial(sys_cl, x0, t);

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
title('Closed-Loop System States Response with Pole Placement');
xlabel('Time (s)');
ylabel('States');
legend('\theta (Arm Angle)', '\alpha (Pendulum Angle)', ' (Arm Angular Velocity)', '(Pendulum Angular Velocity)');
subplot(2, 1, 2);
u = -K_new * x'; % 计算控制输入
plot(t, u);
title('Control Input');
xlabel('Time (s)');
ylabel('u');

% 绘制零极点图
figure;
pzmap(sys_cl);
title('Closed-Loop System Pole-Zero Map with Pole Placement');
