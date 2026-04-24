clc; clear; close all;

%% Ρυθμίσεις
T = 20; dt = 0.01; t = 0:dt:T; N = length(t);
A = [-2.15 0.25; -0.75 -2]; B = [0; 1.5];
u = 1.2*sin(1.5*t) + 0.9*sin(2.8*t) + 0.6*sin(4.2*t) + 0.3*randn(1,N);
omega_fun = @(tk) 0.4 * [sin(tk); cos(0.5*tk)];
Gamma = diag([3.000 0.010 0.010 0.010 0.010 0.010]);

%% Grid σ1/σ2
sigma_vals = linspace(0.001, 9, 360);
mean_errors = zeros(length(sigma_vals));

%% Grid Search
for i = 1:length(sigma_vals)
    for j = 1:length(sigma_vals)
        sigma1 = sigma_vals(i);
        sigma2 = sigma_vals(j);

        x = zeros(2,N); x(:,1) = [1; -1];
        x_hat = zeros(2,N);
        theta_hat = zeros(6,N);
        theta_hat(:,1) = [-2; 0; 0.1; -1.1; -2.3; 1.7];
        ex_norm = zeros(1,N);

        for k = 1:N-1
            omega = omega_fun(t(k));
            dx = A * x(:,k) + B * u(k) + omega;
            x(:,k+1) = x(:,k) + dt * dx;

            phi1 = [x(1,k); x(2,k); u(k); 0; 0; 0];
            phi2 = [0; 0; 0; x(1,k); x(2,k); u(k)];
            y = dx;
            y_hat = [theta_hat(:,k)' * phi1; theta_hat(:,k)' * phi2];
            ey = y - y_hat;

            damping = [sigma1; sigma1; sigma1; sigma2; sigma2; sigma2];
            v = Gamma * (phi1 * ey(1) + phi2 * ey(2) - damping .* theta_hat(:,k));

            % Προβολή περιορισμών
            a11 = theta_hat(1,k); b2 = theta_hat(6,k);
            violates_a11 = ~((-3 < a11) && (a11 < -1));
            violates_b2 = b2 < 1;

            if violates_a11 || violates_b2
                gradg = zeros(6,1);
                if violates_a11, gradg(1) = sign(a11 + 2); end
                if violates_b2, gradg(6) = -1; end

                if gradg' * v <= 0
                    theta_hat(:,k+1) = theta_hat(:,k) + dt * v;
                else
                    proj = v - Gamma * (gradg * gradg') / (gradg' * Gamma * gradg) * v;
                    theta_hat(:,k+1) = theta_hat(:,k) + dt * proj;
                end
            else
                theta_hat(:,k+1) = theta_hat(:,k) + dt * v;
            end

            A_hat = [theta_hat(1,k+1), theta_hat(2,k+1);
                     theta_hat(4,k+1), theta_hat(5,k+1)];
            B_hat = [theta_hat(3,k+1); theta_hat(6,k+1)];
            dx_hat = A_hat * x_hat(:,k) + B_hat * u(k);
            x_hat(:,k+1) = x_hat(:,k) + dt * dx_hat;

            ex_norm(k+1) = norm(x(:,k+1) - x_hat(:,k+1));
        end

        mean_errors(i,j) = mean(ex_norm);
    end
end

%% Top 3 αποτελέσματα
[sorted_errors, idx] = sort(mean_errors(:));
[top_i, top_j] = ind2sub(size(mean_errors), idx(1:3));

fprintf('\n🏆 Κορυφαίοι 3 συνδυασμοί (με Προβολή):\n');
for n = 1:3
    fprintf('σ1 = %.4f, σ2 = %.4f → Μέσο σφάλμα: %.6f\n', ...
        sigma_vals(top_i(n)), sigma_vals(top_j(n)), sorted_errors(n));
end

%% Heatmap
figure;
imagesc(sigma_vals, sigma_vals, mean_errors);
colorbar;
xlabel('\sigma_2'); ylabel('\sigma_1');
title('📊 Μέσο σφάλμα ‖x - \hat{x}‖ με Προβολή');
set(gca, 'YDir', 'normal');
