clc; clear; close all;

%% Πραγματικές παράμετροι
m_true = 0.75;
L_true = 1.25;
c_true = 0.15;
g = 9.81;

A = 4;
omega = 2;
u_func = @(t) A * sin(omega * t);

%% Ρυθμίσεις προσομοίωσης
T_s = 0.1;
T_end = 20;
t = 0:T_s:T_end;
N = length(t);

%% Προσομοίωση του συστήματος (αληθινό q)
x0 = [0; 0];
[t_fine, X_true] = ode45(@(t,x) real_system(t,x,m_true,L_true,c_true,g,u_func), t, x0);

q_clean = X_true(:,1);
u = A * sin(omega * t)';

%% Προσθήκη λευκού θορύβου στο q(t)
sigma_q = 0.05;
rng(1);  % Σταθερός θόρυβος
q_noisy = q_clean + sigma_q * randn(N, 1);

%% Φίλτρο Λ(s) = (s + 1)^2
Lambda_poly = [1 2 1];  % Συντελεστές του (s + 1)^2

D_q    = tf([0 0 1], Lambda_poly);  % q / Λ(s)
D_dq   = tf([0 1 0], Lambda_poly);  % q' / Λ(s)
D_ddq  = tf([1 0 0], Lambda_poly);  % q'' / Λ(s)

%% Εφαρμογή του φίλτρου με lsim
phi1 = lsim(D_ddq, q_noisy, t);  % q'' / Λ(s)
phi2 = lsim(D_dq,  q_noisy, t);  % q' / Λ(s)
phi3 = lsim(D_q,   q_noisy, t);  % q / Λ(s)
Y    = lsim(D_q,   u,       t);  % u / Λ(s)

Phi = [phi1, phi2, phi3];

%% Εκτίμηση παραμέτρων με Ελάχιστα Τετράγωνα
theta_hat = (Phi' * Phi) \ (Phi' * Y);

theta1 = theta_hat(1);  % mL^2
theta2 = theta_hat(2);  % c
theta3 = theta_hat(3);  % mgL

%% Ανάκτηση φυσικών παραμέτρων
L_est = (theta1 * g) / theta3;
m_est = theta1 / L_est^2;
c_est = theta2;

fprintf("Estimated parameters (2β με φίλτρο και θόρυβο):\n");
fprintf("  L = %.4f (true: %.4f)\n", L_est, L_true);
fprintf("  m = %.4f (true: %.4f)\n", m_est, m_true);
fprintf("  c = %.4f (true: %.4f)\n", c_est, c_true);

%% Προσομοίωση με τις εκτιμημένες παραμέτρους
[t_est, X_est] = ode45(@(t,x) real_system(t,x,m_est,L_est,c_est,g,u_func), t, x0);
q_est = X_est(:,1);
e = q_clean - q_est;

%% Γραφικές παραστάσεις
figure;

subplot(3,1,1);
plot(t, q_clean, 'b', 'LineWidth', 1.5); hold on;
plot(t, q_est, 'r--', 'LineWidth', 1.5);
legend('True q(t)', 'Estimated q(t)');
title('2β: True vs Estimated q(t) with Noise + Filter');
ylabel('q(t)');
grid on;

subplot(3,1,2);
plot(t, e, 'k', 'LineWidth', 1.5);
ylabel('Error e(t)');
title('Estimation Error e(t)');
grid on;

subplot(3,1,3);
plot(t, u, 'g', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('u(t)');
title('Input Signal u(t)');
grid on;

%% Συνάρτηση συστήματος
function dxdt = real_system(t, x, m, L, c, g, u_func)
    q = x(1);
    q_dot = x(2);
    u = u_func(t);
    q_ddot = (1 / (m * L^2)) * (u - c * q_dot - m * g * L * q);
    dxdt = [q_dot; q_ddot];
end
