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

A_hat_all = zeros(2,2,N);
B_hat_all = zeros(2,N);
ex_norm = zeros(1,N);
theta_change = zeros(1,N); % ➕ Ταχύτητα μεταβολής παραμέτρων

for k = 1:N-1
    omega = 0.5 * [sin(0.5 * t(k)); cos(0.3 * t(k))];
    dx = A * x(:,k) + B * u(k) + omega;
    x(:,k+1) = x(:,k) + dt * dx;

    phi1 = [x(1,k); x(2,k); u(k); 0; 0; 0];
    phi2 = [0; 0; 0; x(1,k); x(2,k); u(k)];
    y = dx;
    y_hat = [theta_hat(:,k)' * phi1; theta_hat(:,k)' * phi2];
    ey = y - y_hat;

    v = -Gamma * (phi1 * ey(1) + phi2 * ey(2));

    % Προβολή
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

    % Εκτιμήσεις A, B
    A_hat = [theta_hat(1,k+1), theta_hat(2,k+1);
             theta_hat(4,k+1), theta_hat(5,k+1)];
    B_hat = [theta_hat(3,k+1); theta_hat(6,k+1)];
    A_hat_all(:,:,k+1) = A_hat;
    B_hat_all(:,k+1) = B_hat;

    dx_hat = A_hat * x_hat(:,k) + B_hat * u(k);
    x_hat(:,k+1) = x_hat(:,k) + dt * dx_hat;

    ex_norm(k+1) = norm(x(:,k+1) - x_hat(:,k+1));
    theta_change(k+1) = norm(theta_hat(:,k+1) - theta_hat(:,k)); % ➕
end

fprintf('📊 Χωρίς Σ-τροποποίηση - Μέσο σφάλμα: %.6f\n', mean(ex_norm));

%% Γραφήματα

% 1. Καταστάσεις x(t) και εκτιμήσεις
figure;
plot(t, x(1,:), 'b', t, x(2,:), 'g', ...
     t, x_hat(1,:), 'r--', t, x_hat(2,:), 'k--');
legend('x_1', 'x_2', 'x̂_1', 'x̂_2');
title('Καταστάσεις x(t) και εκτιμήσεις \hat{x}(t)');
xlabel('Χρόνος (s)'); ylabel('Τιμή');

% 2. Σφάλμα e_x(t)
figure;
plot(t, x(1,:) - x_hat(1,:), 'r', ...
     t, x(2,:) - x_hat(2,:), 'b');
legend('e_{x1}', 'e_{x2}');
title('Σφάλμα e_x(t) = x(t) - \hat{x}(t)');
xlabel('Χρόνος (s)'); ylabel('Σφάλμα');

% 3. Συνολικό σφάλμα
figure;
plot(t, ex_norm, 'k');
title('‖x(t) − x̂(t)‖ συνολικό σφάλμα');
xlabel('Χρόνος (s)'); ylabel('Νόρμα Σφάλματος');

% 4. Εκτιμήσεις παραμέτρων θ = [a11, a12, b1, a21, a22, b2]
figure;
labels = {'a_{11}', 'a_{12}', 'b_1', 'a_{21}', 'a_{22}', 'b_2'};
for i = 1:6
    subplot(3,2,i);
    plot(t, theta_hat(i,:));
    title(['\theta_', num2str(i), ' (' labels{i} ')']);
    xlabel('Χρόνος (s)');
    ylabel('Τιμή');
end
sgtitle('Εκτιμήσεις παραμέτρων θ(t)');

% 5. Ταχύτητα μεταβολής παραμέτρων
figure;
plot(t, theta_change, 'm');
title('‖θ̂(k+1) − θ̂(k)‖ - Ταχύτητα μεταβολής παραμέτρων');
xlabel('Χρόνος (s)');
ylabel('Μεταβολή');
