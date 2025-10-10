% Linear Step Response
function y_lin = linear_step_response(M,B,x_ss,deltaD,t_vec)
    C = eye(3); D = zeros(3,1);
    sys = ss(M,B,C,D);
    sys_min = minreal(sys,1e-8);
    u = deltaD .* ( (t_vec >= 100) + (t_vec >= 200) + (t_vec >= 300) );
    u = u(:); % ensure column
    y_dev = lsim(sys_min, u, t_vec);
    y_lin = y_dev + repmat(x_ss,length(t_vec),1);
end