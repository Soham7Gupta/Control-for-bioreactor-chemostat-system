%% DNN model for bioreactor states (X, S, P)
clear; clc; close all;

%% Load Data
data = load('bio_reactor_data_new.mat');
t          = data.t_vec(:);     % time vector (N x 1)
x_nl       = data.x_nl;         % nonlinear clean states (N x 3)
x_nl_noisy = data.x_nl_noisy;   % nonlinear noisy states (N x 3)
x_ss       = data.x_ss(:);      % steady-state [X; S; P]
N  = numel(t);
state_names = {'X (biomass)','S (substrate)','P (product)'};

%% Desinging Neural Network
y_data = x_nl_noisy;   % N x 3 (X,S,P)

% Normalizing the inputs and outputs
t_mu    = mean(t);
t_sigma = std(t);
t_scaled = (t - t_mu) / t_sigma;   % N x 1

y_mu    = mean(y_data, 1);         % 1 x 3
y_sigma = std(y_data, [], 1);      % 1 x 3
Y_scaled = (y_data - y_mu) ./ y_sigma;  % N x 3

% ---- Train / validation split ----
Ntrain = round(0.7 * N);          % 70-30 split
idx    = randperm(N);             % shuffle indices
idx_tr = idx(1:Ntrain);
idx_va = idx(Ntrain+1:end);

XTrain = t_scaled(idx_tr, :);     % Ntrain x 1
YTrain = Y_scaled(idx_tr, :);     % Ntrain x 3
XVal   = t_scaled(idx_va, :);     % Nval x 1
YVal   = Y_scaled(idx_va, :);     % Nval x 3

% Defining neural network architecture
numFeatures    = 1;   % only time as input
numResponses   = 3;   % X, S, P
numHiddenUnits = 64;
layers = [
    featureInputLayer(numFeatures, "Name","input")
    fullyConnectedLayer(numHiddenUnits, "Name","fc1")
    reluLayer("Name","relu1")
    fullyConnectedLayer(numHiddenUnits, "Name","fc2")
    reluLayer("Name","relu2")
    fullyConnectedLayer(numHiddenUnits, "Name","fc3")
    reluLayer("Name","relu3")
    fullyConnectedLayer(numResponses, "Name","fc_out")
    regressionLayer("Name","regressionoutput")];

% Training the options
miniBatchSize = 128;

options = trainingOptions('adam', ...
    'MaxEpochs', 500, ...
    'MiniBatchSize', miniBatchSize, ...
    'InitialLearnRate', 1e-3, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{XVal, YVal}, ...
    'ValidationFrequency', 50, ...
    'Plots','training-progress', ...
    'Verbose',true);

%% Network Training
rng(0);
net = trainNetwork(XTrain, YTrain, layers, options);

XAll = t_scaled;                     % N x 1
YPred_scaled = predict(net, XAll);   % N x 3
YPred = YPred_scaled .* y_sigma + y_mu;  % unscale to physical units

%% Calculate R^2 on validation data only
R2_DNN = zeros(3,1);
for i = 1:3
    y_true = x_nl_noisy(idx_va, i);   % noisy plant data (validation)
    y_hat  = YPred(idx_va, i);        % DNN prediction (validation)

    SSE = sum((y_true - y_hat).^2);
    SST = sum((y_true - mean(y_true)).^2);
    R2_DNN(i) = 1 - SSE / SST;
end

disp(' ')
disp('R² (validation): DNN vs Noisy Nonlinear Plant Data')
for i = 1:3
    fprintf('State %d (%s):   R²_DNN_val = %.5f\n', i, state_names{i}, R2_DNN(i));
end

figure('Name','Bioreactor: Nonlinear vs DNN (All States)',...
       'Units','normalized','Position',[0.05 0.05 0.75 0.85]);

for i = 1:3
    subplot(3,1,i);

    % Clean nonlinear
    plot(t, x_nl(:,i), 'g', 'LineWidth', 2.0); hold on;
    % Noisy nonlinear (plant data)
    plot(t, x_nl_noisy(:,i), 'b.', 'MarkerSize', 10);
    % DNN prediction
    plot(t, YPred(:,i), 'm', 'LineWidth', 2.0);
    % Steady state line
    yline(x_ss(i), ':k', 'LineWidth', 1.5);

    xlabel('Time (h)','FontSize',14,'FontWeight','bold');
    ylabel('g/L','FontSize',14,'FontWeight','bold');
    
    title(sprintf('DNN fit: %s   (R squared = %.6f)', ...
        state_names{i}, R2_DNN(i)),...
        'FontSize',15,'FontWeight','bold');

    % Larger legend
    lgd = legend('Clean','Noisy','DNN','SS',...
        'Location','northoutside','Orientation','horizontal');
    set(lgd,'FontSize',12,'FontWeight','bold');

    grid on;
    set(gca,'FontSize',13,'LineWidth',1.2); % Thicker axes lines
    hold off;
end

