clear; clc;

% --- Πραγματικές Παράμετροι ---
m_true = 1.315; b_true = 0.225; k_true = 0.725;

% --- Learning rates ---
gamma_m = 5;
gamma_b = 10;
gamma_k = 5;

% --- Άλλες Ρυθμίσεις ---
a = 10;
tspan = [0 20];
u_type = "sine";

% --- Feedback όροι (μεικτή τοπολογία) ---
theta_m1 = 56;
theta_m2 = 1;

% --- Θόρυβος ---
eta0 = 0.25; f0 = 20;

% --- Αρχικές Συνθήκες ---
x0 = zeros(13,1);
x0(11) = 1.5;   % αρχική μ̂
x0(12) = 0.05;  % αρχική b̂
x0(13) = 0.533; % αρχική k̂

opts = odeset('RelTol',1e-6, 'AbsTol',1e-8);
[t, X] = ode45(@(t,x) lyap_mixed_noise(t,x,m_true,b_true,k_true,...
    a,gamma_m,gamma_b,gamma_k,theta_m1,theta_m2,u_type,eta0,f0), ...
    tspan, x0, opts);

% --- Ανάλυση Αποτελεσμάτων ---
x = X(:,1);
x_hat = X(:,3);
theta_hat = X(:,11:13);
e_x = x - x_hat;

% --- Εκτύπωση Τελικών Εκτιμήσεων ---
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

% --- Συνάρτηση δυναμικής μεικτής δομής με θόρυβο ---
function dx = lyap_mixed_noise(t, x, m, b, k, a, gm, gb, gk, tm1, tm2, u_type, eta0, f0)

    % --- Καταστάσεις ---
    x1 = x(1); x2 = x(2);
    x1_hat = x(3); x2_hat = x(4);
    z0 = x(5); dz0 = x(6);
    z1 = x(7); dz1 = x(8);
    vu = x(9); dvu = x(10);
    m_hat = x(11); b_hat = x(12); k_hat = x(13);

    % --- Είσοδος ---
    switch u_type
        case "sine", u = 2.5 * sin(t);
        otherwise,  u = 0;
    end

    % --- Θόρυβος στην έξοδο ---
    noise = eta0 * sin(2*pi*f0*t);
    x1_noisy = x1 + noise;
    x2_noisy = x2; % μπορείς να προσθέσεις θόρυβο και εδώ αν θέλεις

    % --- Πραγματικό σύστημα ---
    dx1 = x2;
    dx2 = (1/m) * (u - b*x2 - k*x1);

    % --- Σφάλματα με θόρυβο ---
    ex1 = x1_noisy - x1_hat;
    ex2 = x2_noisy - x2_hat;

    % --- Εκτιμητής (μεικτή) ---
    m_hat = max(m_hat, 0.05);
    dx1_hat = x2_hat;
    dx2_hat = (1/m_hat) * (u - b_hat*x2_noisy - k_hat*x1_noisy) + tm1*ex1 + tm2*ex2;

    % --- Φίλτρα ---
    ddz0 = a^2 * x1_noisy - 2*a*dz0 - a^2*z0;
    ddz1 = a^2 * x2_noisy - 2*a*dz1 - a^2*z1;
    ddvu = a^2 * u - 2*a*dvu - a^2*vu;

    % --- Γραμμικός συνδυασμός για m̂ ---
    alpha = 0.8; beta = 0.2;
    e_m = alpha * ex2 + beta * ex1;

    % --- Ενημέρωση παραμέτρων ---
    dm_hat = gm * e_m * vu;
    db_hat = -gb * ex2 * z1;
    dk_hat = -gk * ex1 * z0;

    % --- Περιορισμοί ---
    dm_hat = max(min(dm_hat, 0.5), -0.5);
    db_hat = max(min(db_hat, 0.5), -0.5);
    dk_hat = max(min(dk_hat, 0.5), -0.5);

    m_hat_new = max(min(m_hat + dm_hat, 5), 0.1);
    b_hat_new = max(min(b_hat + db_hat, 5), 0);
    k_hat_new = max(min(k_hat + dk_hat, 5), 0);

    dx = zeros(13,1);
    dx(1:2) = [dx1; dx2];
    dx(3:4) = [dx1_hat; dx2_hat];
    dx(5:6) = [dz0; ddz0];
    dx(7:8) = [dz1; ddz1];
    dx(9:10) = [dvu; ddvu];
    dx(11:13) = [m_hat_new - m_hat; b_hat_new - b_hat; k_hat_new - k_hat];
end
