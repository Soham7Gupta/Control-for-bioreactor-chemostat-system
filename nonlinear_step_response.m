% Non-Linear Step Response
function [t,x] = nonlinear_step_response(f,D0,x_ss,deltaD,t_vec)
    odefun=@(t,x) f(x, D0 + deltaD*(t>=100) + deltaD*(t>=200) + deltaD*(t>=300)); % Used the function f built for fsolve itself
    [t,x] = ode45(odefun,t_vec,x_ss);
end
