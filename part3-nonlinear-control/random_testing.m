clc; clear; close all;

%% 1. Ρυθμίσεις
T = 20; dt = 0.01; t = 0:dt:T; N = length(t);

A = [-2.15 0.25; -0.75 -2];
B = [0; 1.5];

x = zeros(2,N); x_hat = zeros(2,N);
theta1_hat = zeros(3,N); theta2_hat = zeros(3,N);
ex_norm = zeros(1,N);

% ❗ Κακές αρχικές τιμές
theta1_hat(:,1) = [-0.5; -1.0; -1.0];
theta2_hat(:,1) = [-2.5; 0.5; 0.3];

% ❗ Μεγαλύτερο Γ
Gamma = diag([0.6, 0.3, 0.3]);

x(:,1) = [1; -1];
u = 1.2*sin(1.5*t) + 0.9*sin(2.8*t) + 0.6*sin(4.2*t) + 0.3*randn(1,N);

for k = 1:N-1
    dx = A * x(:,k) + B * u(k);
    x(:,k+1) = x(:,k) + dt * dx;
    phi = [x(1,k); x(2,k); u(k)];

    % Εκτιμητής θ1
    y1 = dx(1); y1_hat = theta1_hat(:,k)' * phi; ey1 = y1 - y1_hat;
    v1 = -Gamma * phi * ey1;
    a11_hat = theta1_hat(1,k);
    if (a11_hat > -3) && (a11_hat < -1)
        theta1_hat(:,k+1) = theta1_hat(:,k) + dt * v1;
    else
        gradg = [1; 0; 0] * sign(a11_hat + 2);
        if gradg' * v1 <= 0
            theta1_hat(:,k+1) = theta1_hat(:,k) + dt * v1;
        else
            proj = v1 - Gamma * (gradg * gradg') / (gradg' * Gamma * gradg) * v1;
            theta1_hat(:,k+1) = theta1_hat(:,k) + dt * proj;
        end
    end

    % Εκτιμητής θ2
    y2 = dx(2); y2_hat = theta2_hat(:,k)' * phi; ey2 = y2 - y2_hat;
    v2 = -Gamma * phi * ey2;
    b2_hat = theta2_hat(3,k);
    if b2_hat >= 1
        theta2_hat(:,k+1) = theta2_hat(:,k) + dt * v2;
    else
        gradg = [0; 0; -1];
        if gradg' * v2 <= 0
            theta2_hat(:,k+1) = theta2_hat(:,k) + dt * v2;
        else
            proj = v2 - Gamma * (gradg * gradg') / (gradg' * Gamma * gradg) * v2;
            theta2_hat(:,k+1) = theta2_hat(:,k) + dt * proj;
        end
    end

    % Εκτίμηση x̂
    A_hat = [theta1_hat(1,k+1), theta1_hat(2,k+1);
             theta2_hat(1,k+1), theta2_hat(2,k+1)];
    B_hat = [theta1_hat(3,k+1); theta2_hat(3,k+1)];
    dx_hat = A_hat * x_hat(:,k) + B_hat * u(k);
    x_hat(:,k+1) = x_hat(:,k) + dt * dx_hat;
    ex_norm(k+1) = norm(x(:,k+1) - x_hat(:,k+1));
end

%% Τελικές Εκτιμήσεις
fprintf('Τελική εκτίμηση πίνακα Â:\n');
disp([theta1_hat(1:2,end)'; theta2_hat(1:2,end)']);
fprintf('Τελική εκτίμηση διανύσματος B̂:\n');
disp([theta1_hat(3,end); theta2_hat(3,end)]);

%% Σφάλμα
figure;
plot(t, ex_norm, 'k');
title('Συνολικό σφάλμα ‖x(t) − x̂(t)‖');
xlabel('t'); ylabel('‖e_x(t)‖');
