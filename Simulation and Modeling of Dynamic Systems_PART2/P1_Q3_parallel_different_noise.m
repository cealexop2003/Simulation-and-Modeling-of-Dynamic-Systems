clear; clc;

% --- Πραγματικές Παράμετροι ---
m_true = 1.315; 
b_true = 0.225; 
k_true = 0.725;

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
x0(11) = 0.8;   
x0(12) = 0.2;   
x0(13) = 100;  

opts = odeset('RelTol',1e-6, 'AbsTol',1e-8);

% --- Εύρος πλάτους θορύβου η₀ ---
eta0_values = 0:0.05:0.5;   % Τώρα φτάνει μέχρι 0.5

f0 = 20;                  % Συχνότητα θορύβου (Hz)

% --- Αποθήκευση σφαλμάτων ---
error_m = zeros(size(eta0_values));
error_b = zeros(size(eta0_values));
error_k = zeros(size(eta0_values));

for i = 1:length(eta0_values)
    eta0 = eta0_values(i);

    % --- Λύση χωρίς θόρυβο ---
    [t, X] = ode45(@(t,x) lyap_parallel_dynamics_noise(t,x,m_true,b_true,k_true,...
        a, gamma_m, gamma_b, gamma_k, u_type, eta0, f0), ...
        tspan, x0, opts);

    x = X(:,1);          % Πραγματικό x
    x_hat = X(:,3);      % Εκτιμημένο x_hat
    theta_hat = X(:,11:13); % Εκτιμήσεις παραμέτρων

    % --- Τελικές εκτιμήσεις ---
    m_est = theta_hat(end,1);
    b_est = theta_hat(end,2);
    k_est = theta_hat(end,3);

    % --- Υπολογισμός απόλυτου σφάλματος ---
    error_m(i) = abs(m_est - m_true);
    error_b(i) = abs(b_est - b_true);
    error_k(i) = abs(k_est - k_true);
end

% --- Plot των Σφαλμάτων ---
figure;
plot(eta0_values, error_m, '-o', 'LineWidth', 2); hold on;
plot(eta0_values, error_b, '-s', 'LineWidth', 2);
plot(eta0_values, error_k, '-d', 'LineWidth', 2);
grid on;
xlabel('Πλάτος Θορύβου \eta_0');
ylabel('Σφάλμα Εκτίμησης');
title('Εξέλιξη Σφάλματος Εκτίμησης με το Πλάτος Θορύβου (Παράλληλη)');
legend('|m̂ - m_{true}|', '|b̂ - b_{true}|', '|k̂ - k_{true}|');


function dx = lyap_parallel_dynamics_noise(t, x, m, b, k, a, gamma_m, gamma_b, gamma_k, u_type, eta0, f0)

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

    % --- Θόρυβος μέτρησης ---
    noise = eta0 * sin(2*pi*f0*t);

    % --- Πραγματικό σύστημα με θόρυβο ---
    dx1 = x2;
    dx2 = (1/m) * (u - b*x2 - k*x1);

    x1_noisy = x1 + noise;   % Προσθήκη θορύβου στη μέτρηση

    % --- Εκτιμητής ---
    m_hat = max(m_hat, 0.05);  
    dx1_hat = x2_hat;
    dx2_hat = (1/m_hat) * (u - b_hat*x2_hat - k_hat*x1_hat);

    % --- Σφάλματα ---
    ex1 = x1_noisy - x1_hat;
    ex2 = x2 - x2_hat;

    % --- Φίλτρα ---
    ddz0 = a^2 * x1_hat - 2*a*dz0 - a^2*z0;
    ddz1 = a^2 * x2_hat - 2*a*dz1 - a^2*z1;
    ddvu = a^2 * u      - 2*a*dvu - a^2*vu;

    % --- Εκτίμηση Παραμέτρων ---
    alpha = 0.8;
    beta  = 0.2;
    e_m = alpha * ex2 + beta * ex1;

    dm_hat = -gamma_m * e_m * vu;
    db_hat =  gamma_b * ex1 * z1;
    dk_hat =  gamma_k * ex1 * z0;

    % --- Περιορισμός αλλαγών ---
    dm_hat = max(min(dm_hat, 0.5), -0.5);
    db_hat = max(min(db_hat, 0.5), -0.5);
    dk_hat = max(min(dk_hat, 0.5), -0.5);

    % --- Clamp ---
    m_hat_new = max(min(m_hat + dm_hat, 5), 0.1);
    b_hat_new = max(min(b_hat + db_hat, 5), 0);
    k_hat_new = max(min(k_hat + dk_hat, 5), 0);

    % --- Παράγωγοι Καταστάσεων ---
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
