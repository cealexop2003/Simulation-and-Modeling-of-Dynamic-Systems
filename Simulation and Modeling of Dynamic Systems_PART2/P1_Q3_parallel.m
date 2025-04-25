clear; clc;

% --- Πραγματικές Παράμετροι ---
m_true = 1.315; b_true = 0.225; k_true = 0.725;

% --- Learning rates ---
gamma_m = 10; gamma_b = 10; gamma_k = 10;

% --- Άλλες Ρυθμίσεις ---
a = 15; tspan = [0 20]; u_type = "sine";

% --- Θόρυβος ---
eta0 = 0.25; f0 = 20;   % Αμγ και συχνότητα θορύβου

% --- Αρχικές Συνθήκες ---
x0 = zeros(13,1);
x0(11:13) = [0.8; 0.2; 100];

opts = odeset('RelTol',1e-6, 'AbsTol',1e-8);
[t, X] = ode45(@(t,x) lyap_parallel_noise(t,x,m_true,b_true,k_true,a,...
    gamma_m,gamma_b,gamma_k,u_type,eta0,f0), tspan, x0, opts);

x = X(:,1);
x_hat = X(:,3);
theta_hat = X(:,11:13);
e_x = x - x_hat;

% --- Εκτύπωση αποτελεσμάτων ---
fprintf('\n--- Εκτιμημένες Τελικές Τιμές Παραμέτρων ---\n');
fprintf('m̂ = %.4f, b̂ = %.4f, k̂ = %.4f\n', theta_hat(end,1), theta_hat(end,2), theta_hat(end,3));

% --- Γραφήματα ---
figure;
plot(t, theta_hat(:,1), 'LineWidth', 2); grid on;
xlabel('Χρόνος [s]'); ylabel('m̂(t)');
title('Μεικτή: Εκτίμηση μάζας m̂(t)');

figure;
plot(t, theta_hat(:,2), 'LineWidth', 2); grid on;
xlabel('Χρόνος [s]'); ylabel('b̂(t)');
title('Μεικτή: Εκτίμηση απόσβεσης b̂(t)');

figure;
plot(t, theta_hat(:,3), 'LineWidth', 2); grid on;
xlabel('Χρόνος [s]'); ylabel('k̂(t)');
title('Μεικτή: Εκτίμηση σταθεράς k̂(t)');

figure;
plot(t, x, 'LineWidth', 2); grid on;
xlabel('Χρόνος [s]'); ylabel('x(t)');
title('Κατάσταση x(t) με θόρυβο');

% --- Δυναμική ---
function dx = lyap_parallel_noise(t, x, m, b, k, a, gm, gb, gk, u_type, eta0, f0)

    % --- Καταστάσεις ---
    x1 = x(1); x2 = x(2); x1_hat = x(3); x2_hat = x(4);
    z0 = x(5); dz0 = x(6); z1 = x(7); dz1 = x(8);
    vu = x(9); dvu = x(10); m_hat = x(11); b_hat = x(12); k_hat = x(13);

    % --- Είσοδος ---
    u = 2.5 * sin(t);

    % --- Θόρυβος μέτρησης ---
    noise = eta0 * sin(2*pi*f0*t);
    x1_noisy = x1 + noise;

    % --- Πραγματικό Σύστημα ---
    dx1 = x2;
    dx2 = (1/m) * (u - b*x2 - k*x1);

    % --- Εκτιμητής ---
    m_hat = max(m_hat, 0.05);
    dx1_hat = x2_hat;
    dx2_hat = (1/m_hat) * (u - b_hat*x2_hat - k_hat*x1_hat);

    % --- Σφάλματα (με x1_noisy) ---
    ex1 = x1_noisy - x1_hat;
    ex2 = x2 - x2_hat;

    % --- Φίλτρα ---
    ddz0 = a^2 * x1_hat - 2*a*dz0 - a^2*z0;
    ddz1 = a^2 * x2_hat - 2*a*dz1 - a^2*z1;
    ddvu = a^2 * u      - 2*a*dvu - a^2*vu;

    % --- Ενημέρωση παραμέτρων ---
    alpha = 0.8; beta = 0.2; e_m = alpha*ex2 + beta*ex1;
    dm = -gm * e_m * vu;
    db =  gb * ex1 * z1;
    dk =  gk * ex1 * z0;

    % --- Clamp updates ---
    dm = max(min(dm,0.5), -0.5);
    db = max(min(db,0.5), -0.5);
    dk = max(min(dk,0.5), -0.5);

    % --- Clamp estimates ---
    m_hat_new = max(min(m_hat + dm, 5), 0.1);
    b_hat_new = max(min(b_hat + db, 5), 0);
    k_hat_new = max(min(k_hat + dk, 5), 0);

    dx = zeros(13,1);
    dx(1:2) = [dx1; dx2];
    dx(3:4) = [dx1_hat; dx2_hat];
    dx(5:6) = [dz0; ddz0];
    dx(7:8) = [dz1; ddz1];
    dx(9:10) = [dvu; ddvu];
    dx(11:13) = [m_hat_new - m_hat; b_hat_new - b_hat; k_hat_new - k_hat];
end
