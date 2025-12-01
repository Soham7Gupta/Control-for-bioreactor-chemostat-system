clear; clc; close all;

%% Load identification and simulation data
id  = load('bio_foptd_identification_all.mat');
data = load('bio_reactor_data_new.mat');

t          = data.t_vec(:);
x_nl       = data.x_nl;
x_nl_noisy = data.x_nl_noisy;
x_ss       = data.x_ss(:);
D0         = data.D0;
deltaD     = data.deltaD;

Kp_hat    = id.Kp_hat;
TauP_hat  = id.TauP_hat;
Theta_hat = id.Theta_hat;

n_states    = 3;
state_names = {'X (biomass)','S (substrate)','P (product)'};

%% Common multi-step D(t) profile (same as in optimfunc_bio)
D_vec = D0 * ones(size(t));
D_vec(t >= 100) = D0 - 1*deltaD;
D_vec(t >= 200) = D0 - 2*deltaD;
D_vec(t >= 300) = D0 - 3*deltaD;
u = D_vec - D0;

Ts = t(2) - t(1);

figure(1); clf;

for y_index = 1:n_states

    y_data = x_nl_noisy(:, y_index);
    y_ss   = x_ss(y_index);

    Kp    = Kp_hat(y_index);
    Tau   = TauP_hat(y_index);
    Theta = Theta_hat(y_index);

    % --- FOPTD simulation for this state ---
    Nd = max(0, round(Theta / Ts));
    u_delayed = zeros(size(u));
    if Nd < length(u)
        u_delayed(Nd+1:end) = u(1:end-Nd);
    end

    y_dev = zeros(size(t));
    for k = 1:length(t)-1
        y_dev(k+1) = y_dev(k) + (Ts/Tau)*(-y_dev(k) + Kp*u_delayed(k));
    end

    y_model = y_ss + y_dev;

    % --- Metrics ---
    err  = y_data - y_model;
    SSE  = sum(err.^2);
    ybar = mean(y_data);
    SST  = sum((y_data - ybar).^2);
    R2   = 1 - SSE/SST;

    fprintf('\n===== Validation: D -> %s (state %d) =====\n', ...
            state_names{y_index}, y_index);
    fprintf('  Kp_hat    = %.6g\n', Kp);
    fprintf('  TauP_hat  = %.6g\n', Tau);
    fprintf('  Theta_hat = %.6g\n', Theta);
    fprintf('  SSE       = %.6g\n', SSE);
    fprintf('  R²        = %.6g\n', R2);

    % --- Subplot for this state ---
    subplot(3,1,y_index);
    plot(t, x_nl(:, y_index), 'g', 'LineWidth', 1.0); hold on;
    plot(t, y_data, 'b.', 'MarkerSize', 4);
    plot(t, y_model, 'r--', 'LineWidth', 1.5);
    yline(y_ss, ':k', 'LineWidth', 1);
    xlabel('Time (h)');
    ylabel('State value (g/L)');
    title(sprintf('D \\rightarrow %s  (R^2 = %.3f)', state_names{y_index}, R2));
    legend('Clean','Noisy','FOPTD','SS','Location','best');

    grid on;

end
