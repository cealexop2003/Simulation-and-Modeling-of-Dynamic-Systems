clear; clc;

% --- Πραγματικές Παράμετροι ---
m_true = 1.315; b_true = 0.225; k_true = 0.725;

% --- Learning rates (διαφορετικά) ---
gamma_m = 10;
gamma_b = 10;
gamma_k = 10;

% --- Άλλες Ρυθμίσεις ---
a = 15;
tspan = [0 20];
u_type = "sine";

% --- Αρχικές Συνθήκες ---
x0 = zeros(13,1);
x0(11) = 0.8;   % αρχική μ̂
x0(12) = 0.2;   % αρχική b̂
x0(13) = 100;   % αρχική k̂ (πιο μακριά από 0)


opts = odeset('RelTol',1e-6, 'AbsTol',1e-8);
[t, X] = ode45(@(t,x) lyap_parallel_dynamics(t,x,m_true,b_true,k_true,...
    a, gamma_m, gamma_b, gamma_k, u_type), ...
    tspan, x0, opts);

% --- Ανάλυση Αποτελεσμάτων ---
x = X(:,1); x_hat = X(:,3); theta_hat = X(:,11:13);
e_x = x - x_hat;

% --- Εκτύπωση Τελικών Εκτιμήσεων ---
fprintf('\n--- Εκτιμημένες Τελικές Τιμές Παραμέτρων ---\n');
fprintf('Πραγματικό m = %.4f, Εκτιμημένο m̂ = %.4f\n', m_true, theta_hat(end,1));
fprintf('Πραγματικό b = %.4f, Εκτιμημένο b̂ = %.4f\n', b_true, theta_hat(end,2));
fprintf('Πραγματικό k = %.4f, Εκτιμημένο k̂ = %.4f\n', k_true, theta_hat(end,3));
% --- Γραφήματα Εκτιμήσεων Παραμέτρων ---
figure;
plot(t, theta_hat(:,1), 'LineWidth', 2); grid on;
xlabel('Χρόνος [s]');
ylabel('m̂(t)');
title('Εκτίμηση μάζας m̂(t)');

figure;
plot(t, theta_hat(:,2), 'LineWidth', 2); grid on;
xlabel('Χρόνος [s]');
ylabel('b̂(t)');
title('Εκτίμηση απόσβεσης b̂(t)');

figure;
plot(t, theta_hat(:,3), 'LineWidth', 2); grid on;
xlabel('Χρόνος [s]');
ylabel('k̂(t)');
title('Εκτίμηση σταθεράς ελατηρίου k̂(t)');

% --- Γράφημα κατάστασης x(t) ---
figure;
plot(t, x, 'LineWidth', 2); grid on;
xlabel('Χρόνος [s]');
ylabel('x(t)');
title('Κατάσταση x(t)');
function dx = lyap_parallel_dynamics(t, x, m, b, k, a, gamma_m, gamma_b, gamma_k, u_type)

    % --- Καταστάσεις ---
    x1 = x(1); x2 = x(2);
    x1_hat = x(3); x2_hat = x(4);
    z0 = x(5); dz0 = x(6);
    z1 = x(7); dz1 = x(8);
    vu = x(9); dvu = x(10);
    m_hat = x(11); b_hat = x(12); k_hat = x(13);

    % --- Είσοδος ---
    switch u_type
        case "sine"
            u = 2.5 * sin(t);
        otherwise
            u = 0;
    end

    % --- Πραγματικό σύστημα ---
    dx1 = x2;
    dx2 = (1/m) * (u - b*x2 - k*x1);

    % --- Εκτιμητής ---
    m_hat = max(m_hat, 0.05);  % προστασία από διαίρεση με 0
    dx1_hat = x2_hat;
    dx2_hat = (1/m_hat) * (u - b_hat*x2_hat - k_hat*x1_hat);

    % --- Σφάλματα ---
    ex1 = x1 - x1_hat;
    ex2 = x2 - x2_hat;

    % --- Φίλτρα ---
    ddz0 = a^2 * x1_hat - 2*a*dz0 - a^2*z0;
    ddz1 = a^2 * x2_hat - 2*a*dz1 - a^2*z1;
    ddvu = a^2 * u      - 2*a*dvu - a^2*vu;

    % --- Γραμμικός συνδυασμός ex1 και ex2 για m̂ ---
    alpha = 0.8;
    beta  = 0.2;
    e_m = alpha * ex2 + beta * ex1;

    % --- Εκτίμηση Παραμέτρων με διαγώνιο Γ ---
    dm_hat = -gamma_m * e_m * vu;
    db_hat =  gamma_b * ex1 * z1;
    dk_hat =  gamma_k * ex1 * z0;

    % --- Περιορισμός αλλαγών ---
    dm_hat = max(min(dm_hat, 0.5), -0.5);
    db_hat = max(min(db_hat, 0.5), -0.5);
    dk_hat = max(min(dk_hat, 0.5), -0.5);

    % --- Clamp εκτιμήσεων ---
    m_hat_new = max(min(m_hat + dm_hat, 5), 0.1);
    b_hat_new = max(min(b_hat + db_hat, 5), 0);
    k_hat_new = max(min(k_hat + dk_hat, 5), 0);

    % --- Παράγωγοι καταστάσεων ---
    dx = zeros(13,1);
    dx(1) = dx1;
    dx(2) = dx2;
    dx(3) = dx1_hat;
    dx(4) = dx2_hat;
    dx(5) = dz0;
    dx(6) = ddz0;
    dx(7) = dz1;
    dx(8) = ddz1;
    dx(9) = dvu;
    dx(10) = ddvu;
    dx(11) = m_hat_new - m_hat;
    dx(12) = b_hat_new - b_hat;
    dx(13) = k_hat_new - k_hat;
end
