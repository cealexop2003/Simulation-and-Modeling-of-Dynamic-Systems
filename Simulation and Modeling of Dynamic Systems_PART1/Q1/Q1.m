clc; clear; close all;

% Given parameters
m = 0.75;
L = 1.25;
c = 0.15;
g = 9.81;
A = 4;
omega = 2;

% Time settings
tspan = [0 20]; % Simulation time (0 to 20 sec)
dt = 1e-4; % Time step (less than 10^(-3))
t = 0:dt:20; % Time vector

% Initial conditions [q(0); q_dot(0)]
x0 = [0; 0];

% Solve the system using ode45
[t, X] = ode45(@(t,x) state_equations(t, x, m, L, c, g, A, omega), t, x0);

% Extract states
q = X(:,1);       % q(t)
q_dot = X(:,2);   % q_dot(t)

% Plot states
figure;
subplot(2,1,1);
plot(t, q, 'b', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('q(t)');
title('State q(t)');
grid on;

subplot(2,1,2);
plot(t, q_dot, 'r', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('q̇(t)');
title('State q̇(t)');
grid on;

% Function defining state equations
function dxdt = state_equations(t, x, m, L, c, g, A, omega)
    q = x(1);
    q_dot = x(2);
    u = A * sin(omega * t); % Input function
    
    % State-space representation
    dxdt = zeros(2,1);
    dxdt(1) = q_dot;
    dxdt(2) = (1 / (m * L^2)) * (u - c * q_dot - m * g * L * q);
end
