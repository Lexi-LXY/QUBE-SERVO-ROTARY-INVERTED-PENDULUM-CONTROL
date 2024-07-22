%define parameters
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
Vm=18;
A=[0 0 1 0;0 0 0 1;0 mp^2*Lp^2*Lr*g/(4*JT) 0 0;0 mp*g*Lp*(Jr+mp*Lr^2)/(2*JT) 0 0];
B=[0;0;(Jp+1/4*(mp*Lp^2))/JT;mp*Lp*Lr/(2*JT)];
C=eye(4);
D=zeros(4,1);
%% Check the system's controllability
C_ctrb = [B A*B A*A*B A*A*A*B]
rank_CM = rank(C_ctrb)        % Rank of the controllability matrix.

if (rank_CM == size(C_ctrb,1))
    disp('The controllability matrix has full rank and the system is therefore controllable.');
else
    disp('The controllability matrix does not have full rank and the system is therefore not controllable.');
end
[v,d]=eig(A)
G=ss2tf(A,B,C,D)%transfer function
[NUM DEN]=ss2tf(A,B,C,D)
%torque= k_t*(Vm-km*thetadot)/Rm

%% change the parameters,change the control input T to V_m
B      = (kt/Rm)*B
A(3,3) =A(3,3)-(kt*km/ Rm)*B(3)
A(4,3) =A(4,3)-(kt*km/ Rm)*B(4)


%%
%LQR控制器
Q = diag([10, 10, 10, 10]); % 状态加权矩阵
R = 1; % 输入加权矩阵

K_lqr = lqr(A, B, Q, R);

% 闭环系统矩阵
A_cl = A - B * K_lqr;

% 定义闭环系统
sys_cl = ss(A_cl, B, C, D);


%% 绘制LQR控制器的系统响应图
t = 0:0.01:10; % 时间范围

x0 = [1 * pi / 180; 1 * pi / 180; 0; 0]; % 初始条件

[y, t, x] = initial(sys_cl, x0, t);


%%
state_max = max(x);
state_min = min(x);

% 打印结果
fprintf('Maximum and Minimum values of each state:\n');
for i = 1:size(x, 2)
    fprintf('State %d: Max = %.4f, Min = %.4f\n', i, state_max(i), state_min(i));
end

%%
figure;
subplot(2, 1, 1);
plot(t, x);
title('Closed-Loop System States Response with LQR Control');
xlabel('Time (s)');
ylabel('States');
legend('\theta (Arm Angle)', '\alpha (Pendulum Angle)', ' (Arm Angular Velocity)', '(Pendulum Angular Velocity)');
subplot(2, 1, 2);
u = -K_lqr * x'; % 计算控制输入
plot(t, u);
title('Control Input');
xlabel('Time (s)');
ylabel('u');

%% 绘制零极点图
figure;
pzmap(sys_cl);
title('Closed-Loop System Pole-Zero Map with LQR Control');