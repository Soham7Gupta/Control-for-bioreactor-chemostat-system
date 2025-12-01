
clear; clc; close all;

%% Load step response data from bioreactor simulation
data = load('bio_reactor_data_new.mat');

t        = data.t_vec(:);
x_nl     = data.x_nl;
x_nl_noisy = data.x_nl_noisy;
x_ss     = data.x_ss(:);
deltaD   = data.deltaD;

Ts = t(2) - t(1);

n_states = 3;              % X,S,P
Kp_hat    = zeros(n_states,1);
TauP_hat  = zeros(n_states,1);
Theta_hat = zeros(n_states,1);
SSE       = zeros(n_states,1);
exitflag_all = zeros(n_states,1);

%% Initial guess
Dy_ss_P  = x_nl_noisy(end,3) - x_nl_noisy(1,3);
Kp0      = Dy_ss_P / (deltaD + eps);
Tau0     = (t(end)-t(1))/5;
Theta0   = 0;
x_initial_guess = [Kp0, Tau0, Theta0];

%% Bounds on [Kp, Tau, Theta]
Up_Bound = [  1000,  500, 100 ];
Lo_Bound = [ -1000,  0.1,   0 ];

options = optimset('Display','iter','TolX',1e-10,'TolFun',1e-10);

state_names = {'X (biomass)','S (substrate)','P (product)'};

for y_index = 1:n_states
    fprintf('\n IDENTIFICATION FOR STATE %d: %s \n', ...
            y_index, state_names{y_index});

    % call fmincon with y_index passed into optimfunc_bio
    [opt_x, fval, exitflag] = fmincon( ...
        @(x) optimfunc_bio(x, y_index), ...
        x_initial_guess, [], [], [], [], ...
        Lo_Bound, Up_Bound, [], options);

    Kp_hat(y_index)    = opt_x(1);
    TauP_hat(y_index)  = opt_x(2);
    Theta_hat(y_index) = opt_x(3);
    SSE(y_index)       = fval;
    exitflag_all(y_index) = exitflag;

    fprintf('  Kp_hat    = %.6g\n', opt_x(1));
    fprintf('  TauP_hat  = %.6g\n', opt_x(2));
    fprintf('  Theta_hat = %.6g\n', opt_x(3));
    fprintf('  SSE       = %.6g\n', fval);
end

%% Save all results
save('bio_foptd_identification_all.mat', ...
     'Kp_hat','TauP_hat','Theta_hat','SSE', ...
     'x_ss','t','deltaD','exitflag_all');

disp('FOPTD identification for all three states saved to bio_foptd_identification_all.mat');

