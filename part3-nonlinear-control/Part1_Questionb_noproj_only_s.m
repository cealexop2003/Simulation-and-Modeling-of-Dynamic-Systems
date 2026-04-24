clc; clear; close all;

%% Ρυθμίσεις
T = 20; dt = 0.01; t = 0:dt:T; N = length(t);
A = [-2.15 0.25; -0.75 -2];
B = [0; 1.5];
u = 1.2*sin(1.5*t) + 0.9*sin(2.8*t) + 0.6*sin(4.2*t) + 0.3*randn(1,N);

x = zeros(2,N); x(:,1) = [1; -1];
x_hat = zeros(2,N);
theta_hat = zeros(6,N);
theta_hat(:,1) = [-2; 0; 0.1; -1.1; -2.3; 1.7];
Gamma = diag([3.000 0.010 0.010 0.010 0.010 0.010]);
M = 5; sigma = 0.1;

A_hat_all = zeros(2,2,N);
B_hat_all = zeros(2,N);
ex_norm = zeros(1,N);
theta_norms = zeros(1,N);

for k = 1:N-1
    omega = 0.5 * [sin(t(k)); cos(0.5 * t(k))];
    dx = A * x(:,k) + B * u(k) + omega;
    x(:,k+1) = x(:,k) + dt * dx;

    phi1 = [x(1,k); x(2,k); u(k); 0; 0; 0];
    phi2 = [0; 0; 0; x(1,k); x(2,k); u(k)];
    y = dx;
    y_hat = [theta_hat(:,k)' * phi1; theta_hat(:,k)' * phi2];
    ey = y - y_hat;

    % σ-τροποποίηση
    theta_k = theta_hat(:,k);
    theta_norm = norm(theta_k);
    if theta_norm < M
        sigma_delta = 0;
    elseif theta_norm <= 2*M
        sigma_delta = sigma * (theta_norm/M - 1);
    else
        sigma_delta = sigma;
    end

    v = -Gamma * (phi1 * ey(1) + phi2 * ey(2) + sigma_delta * theta_k);

    % Χωρίς προβολή
    theta_hat(:,k+1) = theta_hat(:,k) + dt * v;

    A_hat = [theta_hat(1,k+1), theta_hat(2,k+1); theta_hat(4,k+1), theta_hat(5,k+1)];
    B_hat = [theta_hat(3,k+1); theta_hat(6,k+1)];
    A_hat_all(:,:,k+1) = A_hat;
    B_hat_all(:,k+1) = B_hat;

    dx_hat = A_hat * x_hat(:,k) + B_hat * u(k);
    x_hat(:,k+1) = x_hat(:,k) + dt * dx_hat;
    ex_norm(k+1) = norm(x(:,k+1) - x_hat(:,k+1));
    theta_norms(k+1) = theta_norm;
end

fprintf('✅ Με Σ-τροποποίηση ΧΩΡΙΣ προβολή - Μέσο σφάλμα: %.6f\n', mean(ex_norm));
