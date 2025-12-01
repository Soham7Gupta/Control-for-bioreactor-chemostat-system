function fun_val = optimfunc_bio(xo, y_index)
    % Load step response data
    data = load('bio_reactor_data_new.mat');

    t          = data.t_vec(:);
    x_nl_noisy = data.x_nl_noisy;
    x_ss       = data.x_ss(:);
    D0         = data.D0;
    deltaD     = data.deltaD;

    % Plant data for selected state
    y_data = x_nl_noisy(:, y_index);
    y_ss   = x_ss(y_index);
    D_vec = D0 * ones(size(t));
    D_vec(t >= 100) = D0 - 1*deltaD;
    D_vec(t >= 200) = D0 - 2*deltaD;
    D_vec(t >= 300) = D0 - 3*deltaD;

    u = D_vec - D0;        % deviation input

    % FOPTD parameters
    Kp    = xo(1);
    Tau   = xo(2);
    Theta = xo(3);

    if Tau <= 0
        fun_val = 1e20;
        return;
    end

    %Discrete first-order
    Ts = t(2) - t(1);
    Nd = max(0, round(Theta / Ts));    % samples of delay

    u_delayed = zeros(size(u));
    if Nd < length(u)
        u_delayed(Nd+1:end) = u(1:end-Nd);
    end

    y_dev = zeros(size(t));
    for k = 1:length(t)-1
        y_dev(k+1) = y_dev(k) + (Ts/Tau)*(-y_dev(k) + Kp*u_delayed(k));
    end

    y_model = y_ss + y_dev;

    % sSSE
    err = y_data - y_model;
    fun_val = sum(err.^2);

end
