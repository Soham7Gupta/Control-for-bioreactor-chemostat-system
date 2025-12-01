function bio_reactor_closed_loop_final()
    clear; clc; close all;
    SetGraphics(); 

    %% Parameters
    Yxs = 0.4;   % g/g
    beta = 0.2;  % h^(-1)
    Pm  = 50;    % g/L
    Ki  = 22;    % g/L
    alpha = 2.2; % g/
    mum = 0.48;  % h^(-1)
    Km  = 1.2;   % g/L
    Sf0 = 20;    % g/L (disturbance variable)
    D0  = 0.202; % h^(-1) (MV dilution rate at SS)

    % Discrete simulation settings (hours)
    sample_t       = 0.05;   % 3 min/sample
    n_steps        = 3000;   % longer run for multi-step analysis
    plot_every     = 10;     % update plots every 10 samples
    pause_time     = 0.05;
    settle_samples = 2;      % Hold everything at steady state for first 2 samples

    % Measurement nois
    rng(0);
    wk = 0.0 * randn(n_steps,1); % No noise

    %% Non-Linear Model
    mu = @(S,P) (mum * (1 - P / Pm) .* S) ./ (Km + S + (S.^2) / Ki);

    f  = @(x,u,Sf_par) [ ...
        -u * x(1) + mu(x(2), x(3)) * x(1);             % dX/dt
         u * (Sf_par - x(2)) - (1 / Yxs) * mu(x(2), x(3)) * x(1); % dS/dt
        -u * x(3) + (alpha * mu(x(2), x(3)) + beta) * x(1)  % dP/dt
    ];
    X_ss = 5.9956 ; 
    S_ss = 5.0109; 
    P_ss = 19.1267;
    x_ss = [X_ss; S_ss; P_ss];
    fprintf('SS at D0=%.6f, Sf=%.3f -> X=%.6f, S=%.6f, P=%.6f\n', ...
        D0, Sf0, X_ss, S_ss, P_ss);

    %% Choose Control Variable
    fprintf('\nSelect controlled variable (CV):\n');
    fprintf('  1 = Biomass X\n');
    fprintf('  2 = Substrate S\n');
    fprintf('  3 = Product P\n');
    CVchoice = input('Choice (1/2/3): ');
    switch CVchoice
        case 1
            y_index  = 1;
            cv_name  = 'X';
            cv_title = 'Biomass X';
            cv_ss    = X_ss;

            % FOPTD parameters for D -> X
            Kp_est    = -37.26;   
            Tau_P_est =  5.601;  
            Tau_D_est =  0.08001;    

        case 2
            y_index  = 2;
            cv_name  = 'S';
            cv_title = 'Substrate S';
            cv_ss    = S_ss;

            % FOPTD parameters for D -> S 
            Kp_est    =  93.16;   
            Tau_P_est =  5.601;    
            Tau_D_est =  0.08001;  

        case 3
            y_index  = 3;
            cv_name  = 'P';
            cv_title = 'Product P';
            cv_ss    = P_ss;

            % FOPTD parameters for D -> P 
            Kp_est    = -148.3;  
            Tau_P_est =  6.781;    
            Tau_D_est =  0.08001;  

        otherwise
            error('Invalid CV choice. Must be 1, 2, or 3.');
    end

    fprintf('\nUsing FOPTD (D -> %s): Kp=%.5g, Tau_P=%.5g h, Tau_D=%.5g h\n', ...
        cv_name, Kp_est, Tau_P_est, Tau_D_est);

    %% Direct-Syntheis Tuning
    Tau_cl = Tau_P_est / 2;   % target closed-loop time constant
    K_c    = (1 / Kp_est) * (Tau_P_est / (Tau_cl + Tau_D_est));  % P gain
    Tau_I  = Tau_P_est;       % integral time
    T_deri = 0.8;             % derivative time

    fprintf('\nSelect mode:\n  1 = Regulatory (multi-step disturbance in Sf)\n  2 = Servo (multi-step setpoint in %s)\n', cv_name);
    Choice1 = input('Choice (1/2): ');
    fprintf('\nSelect controller:\n  1 = P\n  2 = PI\n  3 = PID\n');
    Choice2 = input('Choice (1/2/3): ');

    % Time vector
    Time_vec = (1:n_steps)' * sample_t;

    %% Servo setpoint profile
    sp_seq = cv_ss * ones(n_steps,1);
    step_times  = [300 700 1100 1500 1900 2300 2700];
    sp_offsets  = [+0.5  -0.3  +0.8  -0.6  +0.4  -0.8   0.0];   % offsets (g/L)
    for i = 1:numel(step_times)
        sp_seq(step_times(i):end) = cv_ss + sp_offsets(i);
    end

    %% Regulatory disturbance profile in Sf
    Sf_seq = Sf0 * ones(n_steps,1);
    reg_step_idx   = [500  1500  2500];
    Sf_offsets_abs = [+2   -2    +1];   

    
    for ii = 1:numel(reg_step_idx)
        Sf_seq(reg_step_idx(ii):end) = max(0.1, Sf0 + Sf_offsets_abs(ii));
    end

    % Initializing the values
    y_true   = zeros(n_steps,1);   % true CV
    y_sp     = zeros(n_steps,1);   % setpoint
    D_hist   = zeros(n_steps,1);   % MV trajectory
    Sf_hist  = zeros(n_steps,1);   % disturbance trajectory
    err_hist = zeros(n_steps,1);   % control error

    % ---------------------- Live plotting setup --------------------------
    fig = figure('Name','Closed-loop (Nonlinear Plant) — Live');
    ax1 = subplot(3,1,1); hold(ax1,'on'); grid(ax1,'on');
    h_y  = plot(ax1, Time_vec(1), cv_ss, '-');      % plant output
    h_sp = plot(ax1, Time_vec(1), cv_ss, '-');      % setpoint
    xlabel(ax1,'Time (h)'); ylabel(ax1, sprintf('%s (g/L)', cv_name));
    title(ax1, sprintf('Controlled Variable (%s) vs Setpoint', cv_title));
    legend(ax1, {'Plant output','Setpoint'}, 'Location','best');

    ax2 = subplot(3,1,2); hold(ax2,'on'); grid(ax2,'on');
    h_D  = plot(ax2, Time_vec(1), D0, '-');
    xlabel(ax2,'Time (h)'); ylabel(ax2,'D (h^{-1})');
    title(ax2,'Manipulated Variable (Dilution rate)');

    ax3 = subplot(3,1,3); hold(ax3,'on'); grid(ax3,'on');
    h_Sf = plot(ax3, Time_vec(1), Sf0, '-');
    xlabel(ax3,'Time (h)'); ylabel(ax3,'S_f (g/L)');
    title(ax3,'Disturbance (Feed substrate)');

    xlim(ax1,[0 Time_vec(end)]);
    xlim(ax2,[0 Time_vec(end)]);
    xlim(ax3,[0 Time_vec(end)]);
    drawnow;

    %% CLOSED LOOP SIMULATION
    xk   = x_ss(:);       % state
    Dk   = D0;            % current MV
    D_ss = D0;            % steady-state MV

    % MV saturation
    D_min = -2;  
    D_max = 2;

    % PI/PID accumulators
    Iacc   = 0;        % integral accumulator
    prev_y = cv_ss;    % previous measurement for derivative on PV

    for k = 1:n_steps

        if k <= settle_samples
            % Initial steady-state hold
            y_sp(k) = cv_ss;
            Sf_k    = Sf0;

            [~, xtraj] = ode45(@(t,xx) f(xx, D_ss, Sf_k), [0, sample_t], xk);
            xk = xtraj(end,:).';

            y_true(k)  = xk(y_index);
            D_hist(k)  = D_ss;
            Sf_hist(k) = Sf_k;
            err_hist(k)= 0;

            if (mod(k,plot_every) == 0) || (k <= settle_samples) || (k == n_steps)
                set(h_y,  'XData', Time_vec(1:k), 'YData', y_true(1:k));
                set(h_sp, 'XData', Time_vec(1:k), 'YData', y_sp(1:k));
                set(h_D,  'XData', Time_vec(1:k), 'YData', D_hist(1:k));
                set(h_Sf, 'XData', Time_vec(1:k), 'YData', Sf_hist(1:k));
                drawnow limitrate;
                if pause_time > 0
                    pause(pause_time);
                end
            end
            continue;
        end

        % Mode-dependent setpoint and disturbance
        if Choice1 == 1   % Regulatory
            y_sp(k) = cv_ss;
            Sf_k    = Sf_seq(k);
        else              % Servo
            y_sp(k) = sp_seq(k);
            Sf_k    = Sf0;
        end

        % Plant propagation with current MV
        [~, xtraj] = ode45(@(t,xx) f(xx, Dk, Sf_k), [0, sample_t], xk);
        xk = xtraj(end,:).';

        yk_true    = xk(y_index);           % true CV
        yk         = yk_true + wk(k);       % measured CV
        y_true(k)  = yk_true;
        D_hist(k)  = Dk;
        Sf_hist(k) = Sf_k;

        % Controllers
        err = y_sp(k) - yk;
        err_hist(k) = err;

        dy = 0;   % derivative term (on PV), default 0

        switch Choice2
            case 1  % P controller
                Dk_cmd = D_ss + K_c * err;

            case 2  % PI controller
                Iacc   = Iacc + (sample_t / Tau_I) * err;
                Dk_cmd = D_ss + K_c * (err + Iacc);

            case 3  % PID controller (D on measurement/PV) 
                % derivative on PV: avoids derivative kick on SP steps
                dy     = (prev_y - yk) / sample_t;   % d(-y)/dt
                Iacc   = Iacc + (sample_t / Tau_I) * err;
                Dk_cmd = D_ss + K_c * (err + Iacc + T_deri * dy);
        end

        % update previous measurement for next iteration 
        prev_y = yk;

        %% Saturation + anti-windup 
        Dk_unsat = Dk_cmd;
        Dk = min(max(Dk_cmd, D_min), D_max);

        if Dk ~= Dk_unsat && Choice2 >= 2 && abs(K_c) > 1e-12
            % Back-calculate Iacc to be consistent with saturated Dk
            if Choice2 == 2         % PI
                Iacc = (Dk - D_ss) / K_c - err;
            elseif Choice2 == 3     % PID with D on PV
                Iacc = (Dk - D_ss) / K_c - err - T_deri * dy;
            end
        end
        % ---------------------- Live plot update -------------------------
        if (mod(k,plot_every) == 0) || (k <= settle_samples) || (k == n_steps)
            set(h_y,  'XData', Time_vec(1:k), 'YData', y_true(1:k));
            set(h_sp, 'XData', Time_vec(1:k), 'YData', y_sp(1:k));
            set(h_D,  'XData', Time_vec(1:k), 'YData', D_hist(1:k));
            set(h_Sf, 'XData', Time_vec(1:k), 'YData', Sf_hist(1:k));
            drawnow limitrate;
            if pause_time > 0
                pause(pause_time);
            end
        end
    end
end
function SetGraphics()
    set(0,'DefaultLineLineWidth', 2)
    set(0,'DefaultaxesLineWidth', 2)
    set(0,'DefaultLineMarkerSize',10)
    set(0,'DefaultaxesFontSize', 18)
    set(0,'DefaultTextFontSize', 18)
    set(0,'DefaultaxesFontName', 'arial')
end
