clear
close all
clc



% The following code calls the function nozzle_1d to solve the problem of
% assignment 1, changing the values of alphas for a fixed number of CVs (N)


alpha_p_vect = logspace(-3, -2, 100);
alpha_u_vect = logspace(-3, -2, 100);



m = 1;
N = 50;

toll_u = 1e-6;
toll_p = toll_u;
it_max = 1000;

L=2;       %[m]
rho=1;     %[kg/m^3]
A_in=0.5;  %[m^2]
A_out=0.1; %[m^2]
p_0=10;    %[Pa]
p_out=0;   %[Pa]  NB: pressure is relative
A=@(x) A_in+(A_out- A_in)/L.*x;



errvect_u = zeros(size(alpha_u_vect));
errvect_p = errvect_u;

for i =1:length(alpha_p_vect)
    
    [u, p, x_u, x_p, m_flow] = nozzle_1d(A, N, alpha_p_vect(i), alpha_u_vect(i), toll_u, toll_p, it_max, m);
    
    u_exact= A_out./A(x_u)*sqrt(2*p_0/rho);
    p_exact= p_0.*(1-(A_out./A(x_p)).^2);
    errvect_u(i) = max(abs(u_exact -u));
    errvect_p(i) = max(abs(p_exact - p));
end

% the plot shows that the error decreases with 1/N, as expected

loglog(alpha_p_vect, errvect_p, "-o", alpha_u_vect, errvect_u, '-o', "LineWidth",2)
title("Study of optimal under-relaxation coefficients", "interpreter", "latex")
xlabel("\(\alpha_p\) and \(alpha_u\)", "Interpreter","latex")
ylabel("Error", "interpreter", "latex")
grid on
legend("Pressure error", "Velocity error",  "interpreter", "latex");


%% Now we use the simplex method (Nelder and Mead) to solve this optimization problem
% The number of grid nodes is fixed. When the function nozzle1d_alpha is
% called, this number is printed out
clc
alpha_opt = fminsearch(@nozzle1d_alpha, [0.01, 0.01]);

fprintf("The optimal value for alpha_p is %.5f, while the optimal value for alpha_u is %.5f", alpha_opt(2), alpha_opt(1));

