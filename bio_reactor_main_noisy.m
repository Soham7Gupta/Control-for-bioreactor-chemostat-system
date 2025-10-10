%% Bio- Reactor main file
function bio_reactor_main_noisy()
    clear, clc, close all;
    %% Parameters
    Yxs = 0.4;   % g/g
    beta = 0.2;  % h^(-1)
    Pm = 50;     % g/L
    Ki = 22;     % g/L
    alpha = 2.2; % g/g
    mum = 0.48;  % h^(-1)
    Km = 1.2;    % g/L
    Sf = 20;     % g/L
    D0 = 0.202;  % 
    %% Initial Guess
    x0_guess = [0, 2, 0];
    %% Model Functions
    mu = @(S,p) (mum * (1 - (p / Pm)) * S) / (Km + S + (S.^2) / Ki);
    f = @(x,u) [
        -u * x(1) + mu(x(2), x(3)) * x(1);
        u * (Sf - x(2)) - (1 / Yxs) * mu(x(2), x(3)) * x(1);
        -u * x(3) + (alpha * (mu(x(2), x(3))) + beta) * x(1)
    ];
    opts = optimoptions("fsolve", "Display", "final");
    [x_ss] = fsolve(@(x) f(x, D0), x0_guess, opts);
    fprintf('x_ss = [%.5f, %.5f, %.5f], residual = %.2e', x_ss, norm(f(x_ss, D0)));
    X_ss = x_ss(1); S_ss = x_ss(2); P_ss = x_ss(3);
    %% Derivative Calculations
    A = Km + S_ss + (S_ss.^2 / Ki);
    dA_dS = 1 + 2 * S_ss / Ki;
    % Taking a function G = S / A
    dG_dS = (A - S_ss * dA_dS) / A^2;
    % Partial derivatives of mu
    dmu_dS = mum * (1 - P_ss / Pm) * dG_dS;
    dmu_dP = -(mum / Pm) * (S_ss / A);
    mu_ss = mum * (1 - P_ss / Pm) * S_ss / A;
    %% Jacobian Matrix
    M = zeros(3, 3);
    % f1
    M(1, 1) = -D0 + mu_ss;      % f1_x
    M(1, 2) = X_ss * dmu_dS;    % f1_s
    M(1, 3) = X_ss * dmu_dP;    % f1_p
    % f2
    M(2, 1) = -(1 / Yxs) * mu_ss;          % f2_x
    M(2, 2) = -D0 - (1 / Yxs) * X_ss * dmu_dS;  % f2_s
    M(2, 3) = -(1 / Yxs) * X_ss * dmu_dP;       % f2_p
    % f3
    M(3, 1) = alpha * mu_ss + beta;      % f3_x
    M(3, 2) = alpha * X_ss * dmu_dS;     % f3_s
    M(3, 3) = -D0 + alpha * X_ss * dmu_dP; % f3_p
    % Analytical vector B
    B = [-X_ss;Sf - S_ss;-P_ss];
    fprintf('\nLinearization at (D0=%.6g):\n', D0);
    fprintf('\nJacobian matrix (M) at the steady state is:\n');
    disp(M);
    fprintf('\nInput matrix (B) corresponding to the Jacobian is:\n');
    disp(B);

    %% Transfer Functions
    C = eye(3);     % Creating an identity matrix
    D = zeros(3, 1); % Creating a zero 3x1 matrix
    sys_ss = ss(M, B, C, D);
    sys_tf = tf(sys_ss)
    %% Poles and Zeros of Each Transfer Function
    fprintf('\n---- Poles and Zeros of Individual Transfer Functions ----\n');
    [num_outputs, num_inputs] = size(sys_tf);
    for i = 1:num_outputs
        for j = 1:num_inputs
            % Compute poles and zeros
            p = pole(sys_tf(i, j));
            z = zero(sys_tf(i, j));
            fprintf('\nPoles and zeros of transfer function (%d) are these:\n', i);
            fprintf('Poles:\n');
            disp(p);
            fprintf('Zeros:\n');
            disp(z);
       end
    end
    eigM = eig(M);
    fprintf('Eigenvalues of A:\n'); disp(eigM);
    %% Step response parameters
    deltaD = 0.02 * D0; % Step Change of 2%
    tfinal = 400;
    t_vec = linspace(0, tfinal, 2000);
    % Linear and nonlinear responses
    y_lin = linear_step_response(M, B, x_ss, deltaD, t_vec);
    [t_nl, x_nl] = nonlinear_step_response(f, D0, x_ss, deltaD, t_vec);

    %% Added Noise, Sensor Measurement
    noise_std = 0.01;
    y_lin_noisy = y_lin + noise_std * x_ss .* randn(size(y_lin));
    x_nl_noisy = x_nl + noise_std * x_ss .* randn(size(x_nl));

    %% Plots
    figure("Name", "Linear vs Non-Linear Step Response (with Noise)");
    names = {'X (biomass)', 'S (substrate)', 'P (product)'};
    for i = 1:3
        subplot(3, 1, i);
        plot(t_vec, y_lin_noisy(:, i), 'b-', 'LineWidth', 1.5);
        hold on;
        plot(t_nl, x_nl_noisy(:, i), 'r--', 'LineWidth', 1.2);
        yline(x_ss(i), ':k', 'LineWidth', 1);
        xlabel('Time (h)'); ylabel('g/L');
        title(sprintf('%s: step D = D0 + %.3g (noisy)', names{i}, deltaD));
        legend('Linear prediction', 'Nonlinear simulation', 'x_{ss}', 'Location', 'best');
        grid on; hold off;
    end
end
%%-----------NOTE---------%%
% As we increase deltaD (step change), the model drifts from the nonlinear solution since linearization holds only near steady state.

