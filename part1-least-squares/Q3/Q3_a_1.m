clc; clear; close all;

%% True system parameters
m_true = 0.75;
L_true = 1.25;
c_true = 0.15;
g = 9.81;

A = 4;
omega = 2;
u_func = @(t) A * sin(omega * t);

%% Simulation setup
T_s = 0.1;
t = 0:T_s:20;
N = length(t);
x0 = [0; 0];

%% Simulate true system (clean data)
[~, X] = ode45(@(t,x) real_system(t,x,m_true,L_true,c_true,g,u_func), t, x0);
q = X(:,1);
q_dot = X(:,2);
u = u_func(t)';

%% Add white Gaussian noise
rng(1);                % Fix random seed for reproducibility
sigma = 0.05;          % Standard deviation of noise

q_noisy     = q     + sigma * randn(size(q));
qdot_noisy  = q_dot + sigma * randn(size(q_dot));
u_noisy     = u     + sigma * randn(size(u));

%% Filter setup: Λ(s) = (s + 1)^2
lamda = [1 2 1];
D_q     = tf([0 0 1], lamda);  % 1 / Λ(s)
D_dq    = tf([0 1 0], lamda);  % s / Λ(s)
D_ddq   = tf([1 0 0], lamda);  % s² / Λ(s)

%% Apply filters to noisy signals
q_f     = lsim(D_q, q_noisy, t);       % q / Λ(s)
qdot_f  = lsim(D_q, qdot_noisy, t);    % q̇ / Λ(s)
u_f     = lsim(D_q, u_noisy, t);       % u / Λ(s)
qddot_f = lsim(D_ddq, q_noisy, t);     % q̈ / Λ(s)

%% Least Squares Estimation
Phi = [q_f, qdot_f, u_f];     % Regression matrix
Y = qddot_f;                  % Filtered output
theta_hat = (Phi' * Phi) \ (Phi' * Y);  % LS solution

% Extract intermediate parameters
A21 = theta_hat(1);
A22 = theta_hat(2);
B2  = theta_hat(3);

%% Recover physical parameters
L_est = -g / A21;
mL2 = 1 / B2;
m_est = mL2 / L_est^2;
c_est = -A22 * mL2;

%% Display results
fprintf("3a – Section 2a with noise:\n");
fprintf("  L = %.4f (true: %.4f)\n", L_est, L_true);
fprintf("  m = %.4f (true: %.4f)\n", m_est, m_true);
fprintf("  c = %.4f (true: %.4f)\n", c_est, c_true);

%% Simulate estimated model
[~, X_est] = ode45(@(t,x) real_system(t,x,m_est,L_est,c_est,g,u_func), t, x0);
q_est = X_est(:,1);
e = q - q_est;

%% Plotting
figure;

subplot(3,1,1);
plot(t, q, 'b', 'LineWidth', 1.5); hold on;
plot(t, q_est, 'r--', 'LineWidth', 1.5);
legend('True q(t)', 'Estimated \hat{q}(t)');
ylabel('q(t)');
title('Section 2a with Noise – True vs Estimated q(t)');
grid on;

subplot(3,1,2);
plot(t, e, 'k', 'LineWidth', 1.5);
ylabel('Error e(t)');
title('Estimation Error: e(t) = q(t) - \hat{q}(t)');
grid on;

subplot(3,1,3);
plot(t, u, 'g', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('u(t)');
title('Input Signal u(t)');
grid on;

%% System dynamics
function dxdt = real_system(t, x, m, L, c, g, u_func)
    q = x(1); q_dot = x(2);
    u = u_func(t);
    q_ddot = (1 / (m * L^2)) * (u - c * q_dot - m * g * L * q);
    dxdt = [q_dot; q_ddot];
end
