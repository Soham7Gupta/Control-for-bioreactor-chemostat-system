function results = identify_foptd_bioreactor()
    %% Parameters along with the units
    Yxs   = 0.4;   % g/g
    beta  = 0.2;   % h^(-1)
    Pm    = 50;    % g/L
    Ki    = 22;    % g/L
    alpha = 2.2;   % g/
    mum   = 0.48;  % h^(-1)
    Km    = 1.2;   % g/L
    Sf0   = 20;    % g/L
    D0    = 0.202; % h^(-1)

    %% Model Equations
    mu = @(S,P) (mum * (1 - P / Pm) .* S) ./ (Km + S + (S.^2) / Ki);

    f  = @(x,u,Sf_par) [
        -u * x(1) + mu(x(2), x(3)) * x(1);             % dX/dt
         u * (Sf_par - x(2)) - (1 / Yxs) * mu(x(2), x(3)) * x(1); % dS/dt
        -u * x(3)+ (alpha * mu(x(2), x(3)) + beta) * x(1)  % dP/dt
    ];

    %% Finding the stable steady state
    x0_guess = [0, 2, 0];
    optsFS = optimoptions("fsolve", ...
        "Display","off", ...
        "FunctionTolerance",1e-12, ...
        "StepTolerance",1e-12);

    x_ss = fsolve(@(x) f(x, D0, Sf0), x0_guess, optsFS);
    X_ss = x_ss(1); 
    S_ss = x_ss(2); 
    P_ss = x_ss(3);

    fprintf('\nSteady-State: X=%.4f  S=%.4f  P=%.4f  @ D=%.3f, Sf=%.2f\n', ...
        X_ss, S_ss, P_ss, D0, Sf0);

    %% Linearising the Equations and using Jacobian method
    Aden   = Km + S_ss + (S_ss.^2 / Ki);
    dA_dS  = 1 + 2 * S_ss / Ki;
    dG_dS  = (Aden - S_ss * dA_dS) / Aden^2;
    dmu_dS = mum * (1 - P_ss / Pm) * dG_dS;
    dmu_dP = -(mum / Pm) * (S_ss / Aden);
    mu_ss  = mum * (1 - P_ss / Pm) * S_ss / Aden;

    M = zeros(3,3);
    % f1
    M(1,1) = -D0 + mu_ss;
    M(1,2) = X_ss * dmu_dS;
    M(1,3) = X_ss * dmu_dP;
    % f2
    M(2,1) = -(1 / Yxs) * mu_ss;
    M(2,2) = -D0 - (1 / Yxs) * X_ss * dmu_dS;
    M(2,3) = -(1 / Yxs) * X_ss * dmu_dP;
    % f3
    M(3,1) = alpha * mu_ss + beta;
    M(3,2) = alpha * X_ss * dmu_dS;
    M(3,3) = -D0 + alpha * X_ss * dmu_dP;

    B    = [-X_ss; Sf0 - S_ss; -P_ss];
    C    = eye(3);
    Dmat = zeros(3,1);

    sys_lin  = ss(M,B,C,Dmat);
    sys_lin  = minreal(sys_lin, 1e-7);

    %% Step Simulation
    t = linspace(0, 120, 6000).';      % time vector (h)
    [y_all, ~] = step(sys_lin, t);     % unit step in D
    % y_all: size [6000 x 3], columns = [X, S, P] deviations
    short  = {'X','S','P'};

    % Preallocate for results
    Kp_arr   = zeros(3,1);
    TauP_arr = zeros(3,1);
    TauD_arr = zeros(3,1);
    R2_arr   = zeros(3,1);
    for i = 1:3
        y = y_all(:,i);          % step response (deviation variable)
        Kp = dcgain(sys_lin(i,:));
        % Normalise by Kp so final value ~1
        y_norm = y / Kp;    % should go from 0 -> 1
        % Dead time: first time normalised response reaches 1% of final
        idx_dt = find(y_norm >= 0.01, 1, 'first');
        if isempty(idx_dt)
            Tau_D = 0;
        else
            Tau_D = t(idx_dt);
        end
        % Time constant: time to reach 63.2% of final (after dead time)
        idx_tau = find(y_norm >= 0.632, 1, 'first');
        if isempty(idx_tau) || t(idx_tau) <= Tau_D
            % Fallback: use dominant pole of A matrix
            poles = eig(M);
            [~, idx_dom] = min(abs(real(poles)));
            Tau_P = max(1e-3, -1/real(poles(idx_dom)));
        else
            Tau_P = max(t(idx_tau) - Tau_D, 1e-3);
        end

        % -------- FOPTD step response (deviation) --------
        y_foptd = zeros(size(t));
        for k = 1:length(t)
            if t(k) >= Tau_D
                y_foptd(k) = Kp * (1 - exp(-(t(k) - Tau_D)/Tau_P));
            end
        end

        % -------- R^2 fit --------
        SSE = sum((y - y_foptd).^2);
        SST = sum((y - mean(y)).^2);
        R2  = 1 - SSE/SST;

        % Store in arrays
        Kp_arr(i)   = Kp;
        TauP_arr(i) = Tau_P;
        TauD_arr(i) = Tau_D;
        R2_arr(i)   = R2;
    end
    results = table(Kp_arr, TauP_arr, TauD_arr, R2_arr, ...
        'VariableNames', {'Kp','Tau_P','Tau_D','R2'}, ...
        'RowNames', short);
end
%% Theory for the working mechanism of code to build FOPTD:
%{
Estimating FOPTD parameters (Kp, Tau_D, Tau_P):
   - From the step response of the linearised model, the final steady-state
     change divided by the input step gives the process gain Kp.
   - Then I normalise the response by Kp so it ideally rises from 0 → 1.
   - The "Dead Time" (Tau_D)(delay) is estimated as the first time when the
     normalised response reaches 1% of its final value. This indicates when
     the process first starts reacting to the input change.
   - The time constant (Tau_P) is found by the time gap between Tau_D and
     the moment when the normalised response reaches 63.2% of the final
     value. This represents how fast the system responds after the delay.
   - If the 63.2 point is not clear, we estimate Tau_P using the dominant
     pole of the system (main dynamic mode).
%}