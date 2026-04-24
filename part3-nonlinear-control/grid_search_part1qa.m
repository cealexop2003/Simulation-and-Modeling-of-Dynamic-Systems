clc; clear; close all;

%% Πραγματικό σύστημα
A = [-2.15 0.25; -0.75 -2];
B = [0; 1.5];

%% Χρονικές παραμέτροι
T = 20; dt = 0.01; t = 0:dt:T; N = length(t);
u = 1.2*sin(1.5*t) + 0.9*sin(2.8*t) + 0.6*sin(4.2*t) + 0.3*randn(1,N);

%% Αρχικές συνθήκες
theta0 = [-2; 0; 0.1; -1.1; -2.3; 1.7];  % εντός επιθυμητού συνόλου
x = zeros(2,N); x(:,1) = [1; -1];
x_hat = zeros(2,N);
theta_hat = zeros(6,N); theta_hat(:,1) = theta0;
ex = zeros(1,N);

%% Ευρύτερες τιμές για Γ
gamma_vals = linspace(0.0001, 3, 4);  % 8^6 = 262144 συνδυασμοί
[g1, g2, g3, g4, g5, g6] = ndgrid(gamma_vals, gamma_vals, gamma_vals, gamma_vals, gamma_vals, gamma_vals);
combos = [g1(:), g2(:), g3(:), g4(:), g5(:), g6(:)];
total = size(combos,1);

fprintf("🔍 Συνολικοί συνδυασμοί: %d\n", total);

%% Grid Search
results = [];

for idx = 1:total
    g = combos(idx,:);
    Gamma = diag(g);

    % Reset
    x = zeros(2,N); x(:,1) = [1; -1];
    x_hat = zeros(2,N);
    theta_hat = zeros(6,N); theta_hat(:,1) = theta0;
    ex = zeros(1,N);

    for k = 1:N-1
        dx = A * x(:,k) + B * u(k);
        x(:,k+1) = x(:,k) + dt * dx;
        phi = [x(1,k); x(2,k); u(k)];
        
        y = dx;
        phi_big = [phi', zeros(1,3); zeros(1,3), phi'];
        y_hat = phi_big * theta_hat(:,k);
        e = y - y_hat;
        v = -Gamma * phi_big' * e;

        % Περιορισμοί
        theta_next = theta_hat(:,k) + dt * v;
        theta_next(1) = min(max(theta_next(1), -3), -1);  % a11
        theta_next(6) = max(theta_next(6), 1);            % b2
        theta_hat(:,k+1) = theta_next;

        A_hat = [theta_next(1), theta_next(2); theta_next(4), theta_next(5)];
        B_hat = [theta_next(3); theta_next(6)];
        dx_hat = A_hat * x_hat(:,k) + B_hat * u(k);
        x_hat(:,k+1) = x_hat(:,k) + dt * dx_hat;

        ex(k+1) = norm(x(:,k+1) - x_hat(:,k+1));
    end

    results = [results; mean(ex), g];
end

%% Ταξινόμηση και εμφάνιση
results = sortrows(results, 1);
disp("🏆 Top 5 συνδυασμοί Γ:");
for i = 1:min(5, size(results,1))
    fprintf("\nΣετ #%d\n", i);
    fprintf("Μέσο σφάλμα: %.5f\n", results(i,1));
    fprintf("Γ = diag([%.3f %.3f %.3f %.3f %.3f %.3f])\n", results(i,2:7));
end
