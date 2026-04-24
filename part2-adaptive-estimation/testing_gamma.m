clear; clc;

% Πραγματικές παράμετροι συστήματος
m_true = 1.315;
b_true = 0.225;
k_true = 0.725;
u_type = "sine";
tspan = [0 20];

% Σταθερές επιλογές (καλές τιμές από προηγούμενο grid)
a = 15;
omega_n = a;
zeta = 1;
theta0 = [0.5; 0.2; 0.7];

% Λίστα για πολλά Gamma (20+ διαφορετικά combos)
gamma_vals = [
    5   5   5;
    10  10  10;
    15  10  10;
    10  15  10;
    10  10  15;
    20  10  10;
    10  20  10;
    10  10  20;
    15  15  10;
    10  15  15;
    15  10  15;
    20  20  20;
    25  10  10;
    10  25  10;
    10  10  25;
    15  20  10;
    10  15  20;
    20  15  15;
    15  15  15;
    30  30  30
];

% Αποθήκευση αποτελεσμάτων
results = [];

% Loop σε όλα τα Γάμμα
for i = 1:size(gamma_vals,1)
    Gamma = diag(gamma_vals(i,:));
    
    % Αρχικές συνθήκες
    x0 = zeros(13,1);
    x0(11:13) = theta0;

    % Προσομοίωση
    [t, X] = ode45(@(t,x) dynamics(t, x, m_true, b_true, k_true, Gamma, omega_n, zeta, u_type), tspan, x0);
    theta_hat = X(:,11:13);
    theta_final = theta_hat(end,:);

    % Σφάλμα
    err = abs(theta_final - [m_true, b_true, k_true]);
    total_err = sum(err);

    % Καταγραφή
    results(end+1,:) = [gamma_vals(i,:), theta_final, total_err]; %#ok<AGROW>
end

% Ταξινόμηση
results_sorted = sortrows(results, 7);

% Εμφάνιση
fprintf('\n%-18s   %-18s   %-10s\n', 'Gamma = diag(...)', 'θ̂(t_end)', 'TotalErr');
for i = 1:size(results_sorted,1)
    row = results_sorted(i,:);
    fprintf('[%4.0f %4.0f %4.0f]   [%5.3f %5.3f %5.3f]   %.4f\n', ...
        row(1:3), row(4:6), row(7));
end

% Καλύτερο run
best = results_sorted(1,:);
best_gamma = diag(best(1:3));
x0 = zeros(13,1); x0(11:13) = theta0;

[t, X] = ode45(@(t,x) dynamics(t, x, m_true, b_true, k_true, best_gamma, omega_n, zeta, u_type), tspan, x0);
theta_hat = X(:,11:13);

figure;
plot(t, theta_hat, 'LineWidth', 1.5); grid on;
xlabel('Time [s]'); ylabel('Estimated Parameters');
legend('\hat{m}(t)', '\hat{b}(t)', '\hat{k}(t)');
title(sprintf('Best Gamma = diag([%.0f %.0f %.0f])', best(1:3)));

function dx = dynamics(t, x, m, b, k, Gamma, omega_n, zeta, u_type)

% Καταστάσεις
x1 = x(1); x2 = x(2);
z = x(3); dz = x(4);
z1 = x(5); dz1 = x(6);
z2 = x(7); dz2 = x(8);
vu = x(9); dvu = x(10);
theta_hat = x(11:13);

% Είσοδος
switch u_type
    case "const"
        u = 2.5;
    case "sine"
        u = 2.5 * sin(t);
    otherwise
        u = 0;
end

% Δυναμική πραγματικού συστήματος
ddx = (1/m)*(u - b*x2 - k*x1);

% Φίλτρα 2ης τάξης: Λ(s) = (s + a)^2
ddz  = omega_n^2*x1  - 2*zeta*omega_n*dz  - omega_n^2*z;
ddz1 = omega_n^2*x2  - 2*zeta*omega_n*dz1 - omega_n^2*z1;
ddz2 = omega_n^2*ddx - 2*zeta*omega_n*dz2 - omega_n^2*z2;
ddvu = omega_n^2*u   - 2*zeta*omega_n*dvu - omega_n^2*vu;

% Gradient update
phi = [z2; z1; z];
v = vu;
e = phi' * theta_hat - v;
dtheta_hat = -Gamma * phi * e;

% Παράγωγοι καταστάσεων
dx = zeros(13,1);
dx(1) = x2;
dx(2) = ddx;
dx(3) = dz;     dx(4) = ddz;
dx(5) = dz1;    dx(6) = ddz1;
dx(7) = dz2;    dx(8) = ddz2;
dx(9) = dvu;    dx(10) = ddvu;
dx(11:13) = dtheta_hat;

end
