clc; clear; close all;

%% 1. Ρυθμίσεις
T = 20; dt = 0.01; t = 0:dt:T; N = length(t);

%% 2. Είσοδος
u = 1.2*sin(1.5*t) + 0.9*sin(2.8*t) + 0.6*sin(4.2*t) + 0.3*randn(1,N);

%% 3. Τιμές πλέγματος (grid) για θ0 και Γ
theta_vals = linspace(-2, 1, 2);        % 3 τιμές για κάθε θ0 συνιστώσα
gamma_vals = linspace(0.01, 1, 3);      % 3 τιμές για κάθε γ (συνολικά 3^6 = 729 συνδυασμοί)

results = []; % για αποθήκευση αποτελεσμάτων
combo_id = 0;

%% 4. Ορισμός πραγματικού συστήματος
A = [-2.15 0.25; -0.75 -2];
B = [0; 1.5];

%% 5. Όλοι οι συνδυασμοί
for a11_0 = theta_vals
for a12_0 = theta_vals
for b1_0  = theta_vals
for a21_0 = theta_vals
for a22_0 = theta_vals
for b2_0  = theta_vals

    theta0 = [a11_0; a12_0; b1_0; a21_0; a22_0; b2_0];

    % Παράλειψη εάν είμαστε πολύ κοντά στη σωστή λύση
    if all(abs(theta0 - [-2.15; 0.25; 0; -0.75; -2; 1.5]) < 0.2)
        continue;
    end

    for g1 = gamma_vals
    for g2 = gamma_vals
    for g3 = gamma_vals
    for g4 = gamma_vals
    for g5 = gamma_vals
    for g6 = gamma_vals

        combo_id = combo_id + 1;
        Gamma = diag([g1, g2, g3, g4, g5, g6]);

        %% Αρχικοποιήσεις
        x = zeros(2,N); x(:,1) = [1; -1];
        x_hat = zeros(2,N);
        theta_hat = zeros(6,N); theta_hat(:,1) = theta0;
        ex_norm = zeros(1,N);

        for k = 1:N-1
            dx = A * x(:,k) + B * u(k);
            x(:,k+1) = x(:,k) + dt * dx;

            % Φ διανύσματα
            phi1 = [x(1,k); x(2,k); u(k); 0; 0; 0];
            phi2 = [0; 0; 0; x(1,k); x(2,k); u(k)];
            phi = phi1 + phi2;

            % Έξοδοι
            y = dx;
            y_hat = [theta_hat(:,k)' * phi1; theta_hat(:,k)' * phi2];
            ey = y - y_hat;

            % Νόμος προσαρμογής
            v = -Gamma * (phi1 * ey(1) + phi2 * ey(2));

            % Περιορισμοί
            a11 = theta_hat(1,k); b2 = theta_hat(6,k);
            violates_a11 = ~((-3 < a11) && (a11 < -1));
            violates_b2 = b2 < 1;

            % Προβολή αν χρειάζεται
            if violates_a11 || violates_b2
                gradg = zeros(6,1);
                if violates_a11
                    gradg(1) = sign(a11 + 2);
                end
                if violates_b2
                    gradg(6) = -1;
                end
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

            % Εκτίμηση κατάστασης
            dx_hat = A_hat * x_hat(:,k) + B_hat * u(k);
            x_hat(:,k+1) = x_hat(:,k) + dt * dx_hat;

            % Σφάλμα
            ex_norm(k+1) = norm(x(:,k+1) - x_hat(:,k+1));
        end

        % Μέσο σφάλμα
        mean_err = mean(ex_norm);
        results = [results; mean_err, theta0', g1, g2, g3, g4, g5, g6];

    end; end; end; end; end; end
end; end; end; end; end; end

%% 6. Αποτελέσματα
results = sortrows(results, 1);
disp("🏁 Top 5 συνδυασμοί (θ₀ και Γ):");
for i = 1:min(5, size(results,1))
    fprintf("\nΣετ #%d:\n", i);
    fprintf("Μέσο σφάλμα: %.5f\n", results(i,1));
    fprintf("θ₀ = [%g, %g, %g, %g, %g, %g]\n", results(i,2:7));
    fprintf("Γ = diag([%.3f, %.3f, %.3f, %.3f, %.3f, %.3f])\n", results(i,8:13));
end
