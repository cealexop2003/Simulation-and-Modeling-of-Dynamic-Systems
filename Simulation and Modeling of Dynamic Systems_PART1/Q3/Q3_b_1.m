clc; clear; close all;

% True parameters
m_true = 0.75; L_true = 1.25; c_true = 0.15; g = 9.81;
A = 4; omega = 2;
u_func = @(t) A * sin(omega * t);
lamda = [1 2 1];  % Filter

% Transfer functions
D_q = tf([0 0 1], lamda);
D_dq = tf([0 1 0], lamda);
D_ddq = tf([1 0 0], lamda);

% Sampling periods to test
Ts_values = 0.01:0.01:0.3;
nTs = length(Ts_values);
errors = zeros(nTs, 3);  % [L, m, c]

for i = 1:nTs
    T_s = Ts_values(i);
    t = 0:T_s:20;
    x0 = [0; 0];

    % Simulate system
    [~, X] = ode45(@(t,x) real_system(t,x,m_true,L_true,c_true,g,u_func), t, x0);
    q = X(:,1); q_dot = X(:,2);
    u = u_func(t)';

    % Filtered signals
    q_f = lsim(D_q, q, t);
    qdot_f = lsim(D_q, q_dot, t);
    u_f = lsim(D_q, u, t);
    qddot_f = lsim(D_ddq, q, t);

    % Least squares
    Phi = [q_f, qdot_f, u_f];
    Y = qddot_f;
    theta = (Phi' * Phi) \ (Phi' * Y);

    A21 = theta(1); A22 = theta(2); B2 = theta(3);

    % Recover parameters
    L_est = -g / A21;
    mL2 = 1 / B2;
    m_est = mL2 / L_est^2;
    c_est = -A22 * mL2;

    % Relative errors (%)
    errors(i,:) = abs([L_est, m_est, c_est] - [L_true, m_true, c_true]) ./ [L_true, m_true, c_true] * 100;
end

% Plotting
params = {'L', 'm', 'c'};
for j = 1:3
    figure;
    plot(Ts_values, errors(:,j), 'b-o', 'LineWidth', 1.5);
    xlabel('Sampling Period T_s (s)');
    ylabel(['Relative Error (%) in ', params{j}]);
    title(['Section 3b – 2a: Error vs T_s for ', params{j}]);
    grid on;
end

% System function
function dxdt = real_system(t, x, m, L, c, g, u_func)
    q = x(1); q_dot = x(2); u = u_func(t);
    q_ddot = (1 / (m * L^2)) * (u - c * q_dot - m * g * L * q);
    dxdt = [q_dot; q_ddot];
end


