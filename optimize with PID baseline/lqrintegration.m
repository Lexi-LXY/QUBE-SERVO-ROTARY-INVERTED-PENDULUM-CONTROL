%define parameters
rng(22);
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


%% change the parameters,change the control input T to V_m
B      = (kt/Rm)*B
A(3,3) =A(3,3)-(kt*km/ Rm)*B(3)
A(4,3) =A(4,3)-(kt*km/ Rm)*B(4)




%%  PSO 优化设置
%% 时间步长和仿真时间
% PSO optimization settings
Ts = 0.01; % Time step
T = 0:Ts:10; % Simulation time
x0 = [1 * pi / 180; 1 * pi / 180; 0; 0]; % Initial condition
target_integral_abs = 0.0473; % Target control input integral absolute value

% Define objective function
objectiveFunction = @(params) objectiveFunctionPSO(params, A, B, Ts, target_integral_abs, T, x0);

% Initial guess for Q and R
initialParams = [10, 10, 10, 10, 1];

% PSO optimization
options = optimoptions('particleswarm', 'SwarmSize', 30, 'MaxIterations', 300, 'Display', 'iter');
[params_opt, fval] = particleswarm(objectiveFunction, 5, [0, 0, 0, 0, 0], [100, 100, 100, 100, 100], options);

% Extract optimized Q and R
Q_opt = diag(params_opt(1:4));
R_opt = params_opt(5);

disp('Optimized Q:');
disp(Q_opt);
disp('Optimized R:');
disp(R_opt);
% Calculate optimized LQR gain
[K_opt, ~, ~] = lqr(A, B, Q_opt, R_opt);

% Define optimized closed-loop system
A_cl_opt = A - B * K_opt;
sys_cl_opt = ss(A_cl_opt, B, C, D);

% Simulate optimized LQR controller
[y_opt, t_opt, x_opt] = initial(sys_cl_opt, x0, T);
u_opt = -K_opt * x_opt';

% Calculate optimized control input integral
integral_opt = sum(abs(u_opt)) * Ts;
fprintf('Optimized LQR Control Input Integral: %.4f\n', integral_opt);

% Plot system response
figure;
subplot(2, 1, 1);
plot(t_opt, x_opt);
title('Closed-Loop System States Response with Optimized LQR Control');
xlabel('Time (s)');
ylabel('States');
legend('\theta (Arm Angle)', '\alpha (Pendulum Angle)', ' (Arm Angular Velocity)', '(Pendulum Angular Velocity)');

subplot(2, 1, 2);
plot(t_opt, u_opt);
title('Control Input');
xlabel('Time (s)');
ylabel('u');

% Plot pole-zero map
figure;
pzmap(sys_cl_opt);
title('Closed-Loop System Pole-Zero Map with Optimized LQR Control');

% Print maximum and minimum values of each state
state_max = max(x_opt);
state_min = min(x_opt);

fprintf('Maximum and Minimum values of each state:\n');
for i = 1:size(x_opt, 2)
    fprintf('State %d: Max = %.4f, Min = %.4f\n', i, state_max(i), state_min(i));
end

% Objective function for PSO
function J = objectiveFunctionPSO(params, A, B, Ts, target_integral_abs, T, x0)
    % Extract Q and R matrix parameters
    Q = diag(params(1:4));
    R = params(5);

    % Ensure R is positive definite
    if R <= 0
        J = Inf;
        return;
    end
    
    % Calculate LQR gain
    try
        [K, ~, ~] = lqr(A, B, Q, R);
    catch
        J = Inf;
        return;
    end
    
    % Closed-loop system matrix
    A_cl = A - B * K;

    % Define closed-loop system
    sys_cl = ss(A_cl, B, eye(4), zeros(4, 1));

    % Simulate closed-loop system
    [~, ~, x] = initial(sys_cl, x0, T);
    u = -K * x';

    % Calculate control input integral
    integral_u = sum(abs(u)) * Ts;

    % Objective function: absolute difference between control input integral and target value
    J = abs(integral_u - target_integral_abs);
end