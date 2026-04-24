% lab2_gradient_estimator.m
% Εκτίμηση παραμέτρων m, b, k με μέθοδο κλίσης και φίλτρο Λ(s) = (s+a)^2

clear; clc;

%% --- ΠΑΡΑΜΕΤΡΟΙ ΣΥΣΤΗΜΑΤΟΣ ---
m_true = 1.315;
b_true = 0.225;
k_true = 0.725;

Gamma = diag([25, 10, 10]);  % Learning rate (τονική προσαρμογή για m)

a = 15;                      % Παράμετρος φίλτρου
omega_n = a; zeta = 1;       % Κλασική μορφή για φίλτρο 2ης τάξης

u_type = "sine";             % ή "const"
tspan = [0 20];

%% --- ΑΡΧΙΚΕΣ ΣΥΝΘΗΚΕΣ ---
x0 = zeros(13,1);
x0(11:13) = [0.5; 0.2; 0.7];

%% --- ΠΡΟΣΟΜΟΙΩΣΗ ---
[t, X] = ode45(@(t,x) dynamics(t, x, m_true, b_true, k_true, Gamma, omega_n, zeta, u_type), tspan, x0);

% Ανάκτηση μεταβλητών
x_true = X(:,1);
theta_hat = X(:,11:13);
m_hat = theta_hat(:,1);
b_hat = theta_hat(:,2);
k_hat = theta_hat(:,3);

%% --- ΓΡΑΦΗΜΑΤΑ ---

% Εκτίμηση m̂(t)
figure;
plot(t, m_hat, 'LineWidth', 1.5); grid on;
xlabel('Time [s]');
ylabel('m̂(t)');
title('Εκτίμηση της μάζας m̂(t)');

% Εκτίμηση b̂(t)
figure;
plot(t, b_hat, 'LineWidth', 1.5); grid on;
xlabel('Time [s]');
ylabel('b̂(t)');
title('Εκτίμηση της απόσβεσης b̂(t)');

% Εκτίμηση k̂(t)
figure;
plot(t, k_hat, 'LineWidth', 1.5); grid on;
xlabel('Time [s]');
ylabel('k̂(t)');
title('Εκτίμηση της ελαστικότητας k̂(t)');

% Θέση x(t)
figure;
plot(t, x_true, 'LineWidth', 1.5); grid on;
xlabel('Time [s]');
ylabel('x(t)');
title('Θέση του Συστήματος x(t)');

%% --- ΕΚΤΥΠΩΣΗ ΤΕΛΙΚΩΝ ΤΙΜΩΝ ---
fprintf('\n--- Εκτιμημένες Τελικές Τιμές Παραμέτρων ---\n');
fprintf('Πραγματικό m = %.4f, Εκτιμημένο m̂ = %.4f\n', m_true, m_hat(end));
fprintf('Πραγματικό b = %.4f, Εκτιμημένο b̂ = %.4f\n', b_true, b_hat(end));
fprintf('Πραγματικό k = %.4f, Εκτιμημένο k̂ = %.4f\n', k_true, k_hat(end));

fprintf('\n--- Απόλυτο Σφάλμα Εκτίμησης ---\n');
fprintf('|m̂ - m| = %.4e\n', abs(m_hat(end) - m_true));
fprintf('|b̂ - b| = %.4e\n', abs(b_hat(end) - b_true));
fprintf('|k̂ - k| = %.4e\n', abs(k_hat(end) - k_true));

%% --- ΔΥΝΑΜΙΚΗ ΣΥΣΤΗΜΑΤΟΣ ---
function dx = dynamics(t, x, m, b, k, Gamma, omega_n, zeta, u_type)

% --- Καταστάσεις ---
x1 = x(1);        % x(t)
x2 = x(2);        % dx/dt
z = x(3); dz = x(4);
z1 = x(5); dz1 = x(6);
z2 = x(7); dz2 = x(8);
vu = x(9); dvu = x(10);
theta_hat = x(11:13); % [m̂, b̂, k̂]

% --- Είσοδος ---
switch u_type
    case "const"
        u = 2.5;
    case "sine"
        u = 2.5 * sin(t);
    otherwise
        u = 0;
end

% --- Πραγματική δυναμική συστήματος ---
ddx = (1/m)*(u - b*x2 - k*x1);

% --- Φίλτρα 2ης τάξης: Λ(s) = (s + a)^2 ---
ddz  = omega_n^2*x1  - 2*zeta*omega_n*dz  - omega_n^2*z;
ddz1 = omega_n^2*x2  - 2*zeta*omega_n*dz1 - omega_n^2*z1;
ddz2 = omega_n^2*ddx - 2*zeta*omega_n*dz2 - omega_n^2*z2;
ddvu = omega_n^2*u   - 2*zeta*omega_n*dvu - omega_n^2*vu;

% --- Gradient Update ---
phi = [z2; z1; z];
v = vu;
e = phi' * theta_hat - v;
dtheta_hat = -Gamma * phi * e;

% --- Παράγωγοι καταστάσεων ---
dx = zeros(13,1);
dx(1) = x2;
dx(2) = ddx;
dx(3) = dz;     dx(4) = ddz;
dx(5) = dz1;    dx(6) = ddz1;
dx(7) = dz2;    dx(8) = ddz2;
dx(9) = dvu;    dx(10) = ddvu;
dx(11:13) = dtheta_hat;

end
