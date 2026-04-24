clc; clear; close all;

% Πραγματικοί πίνακες
A = [-2.15 0.25; -0.75 -2];
B = [0; 1.5];

% Παράμετροι
T = 20; dt = 0.01; t = 0:dt:T; N = length(t);
u = 1.2*sin(1.5*t) + 0.9*sin(2.8*t) + 0.6*sin(4.2*t) + 0.3*randn(1,N);
num_tests = 500;
real_theta = [-2.15 0.25 0 -0.75 -2 1.5];
gamma_vals = [0.05, 0.1, 0.2, 0.3];

results = [];

for test = 1:num_tests
    % --- Τυχαίες αρχικές τιμές θ1, θ2
    theta1_0 = [-3 + 5*rand(), -1 + 2*rand(), -1 + 2*rand()];
    theta2_0 = [-1.5 + 1.2*rand(), -2.8 + 1.6*rand(), 0.4 + 1.2*rand()];
    theta_0 = [theta1_0, theta2_0];

    % Απόρριψη αν κοντά στις πραγματικές
    if any(abs(theta_0 - real_theta) < 0.2)
        continue
    end

    % --- Τυχαίο Γ
    g = datasample(gamma_vals, 3);
    Gamma = diag(g);

    % --- Εκτελούμε εξομοίωση
    x = zeros(2,N); x_hat = zeros(2,N);
    theta1_hat = zeros(3,N); theta2_hat = zeros(3,N);
    ex_norm = zeros(1,N);

    x(:,1) = [1; -1];
    theta1_hat(:,1) = theta1_0';
    theta2_hat(:,1) = theta2_0';

    for k = 1:N-1
        dx = A * x(:,k) + B * u(k);
        x(:,k+1) = x(:,k) + dt * dx;
        phi = [x(1,k); x(2,k); u(k)];

        % Θ1
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

        % Θ2
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

        % x̂
        A_hat = [theta1_hat(1,k+1), theta1_hat(2,k+1);
                 theta2_hat(1,k+1), theta2_hat(2,k+1)];
        B_hat = [theta1_hat(3,k+1); theta2_hat(3,k+1)];
        dx_hat = A_hat * x_hat(:,k) + B_hat * u(k);
        x_hat(:,k+1) = x_hat(:,k) + dt * dx_hat;

        ex_norm(k+1) = norm(x(:,k+1) - x_hat(:,k+1));
    end

    % Καταγραφή σφάλματος και ρυθμίσεων
    mean_err = mean(ex_norm);
    results = [results; mean_err, theta_0, g];
end

% Ταξινόμηση
results = sortrows(results, 1);

% Εμφάνιση Top 5
disp("🏆 Top 5 αποτελέσματα:");
for i = 1:min(5, size(results,1))
    fprintf("\nΣετ #%d\n", i);
    fprintf("Μέσο σφάλμα: %.5f\n", results(i,1));
    fprintf("θ1(0) = [%g, %g, %g]\n", results(i,2:4));
    fprintf("θ2(0) = [%g, %g, %g]\n", results(i,5:7));
    fprintf("Γ = diag([%.3f, %.3f, %.3f])\n", results(i,8:10));
end
