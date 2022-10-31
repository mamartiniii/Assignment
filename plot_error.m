clear
close all
clc

% The following code calls the function nozzle_1d to solve the problem of
% assignment 1, changing the number of nodes (N, hence changing the
% dimension of each CV), and the under-relaxation coefficients.


m = 1;
N = 21;
alpha_p = 0.01;%0.00934; %0.00132;  %these are the values that minimize the error for N = 1000
alpha_u = 0.01;%0.03460; %0.08068;
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

NN = 10:10:80;
hh = 1./NN;
errvect_u = zeros(size(hh));
errvect_p = zeros(size(hh));
m_flow_vect = errvect_p;


for i =1:length(NN)
    N = NN(i);
    [u, p, x_u, x_p, m_flow] = nozzle_1d(A, N, alpha_p, alpha_u, toll_u, toll_p, it_max, m);
    m_flow_vect(i) = m_flow;
    u_exact= A_out./A(x_u)*sqrt(2*p_0/rho);
    p_exact= p_0.*(1-(A_out./A(x_p)).^2);
    errvect_u(i) = max(abs(u_exact -u));
    errvect_p(i) = max(abs(p_exact - p));
end

% the plot shows that the error decreases with 1/N, as expected

loglog(NN, errvect_p, "-o", NN, errvect_u, '-o', NN, 1./NN, "--", NN, 1./(NN.^2), "--", "LineWidth",3)
title("Empirical estimation of order of accuracy of the SIMPLE algorithm", "interpreter", "latex", FontSize=25)
xlabel("Number of nodes \(N\)", "Interpreter","latex", FontSize=20)
ylabel("Error", "interpreter", "latex", fontsize=20)
grid on
legend("\(err_p\)", "\(err_u\)", "\(1/N\)", "\(1/N^2\)",  "interpreter", "latex", fontsize=20);
set(gca,'FontSize',15)
fprintf("The SIMPLE method is a 1st order method (due to the QUICK algorithm, which uses the upwind approximation)\n")


%% Here we compare the exact value of mass flow rate with the one obtained numerically varying N 
 
m_exact =  rho*A_in*A_out/A_in*sqrt(2*p_0/rho);

figure
plot(NN, m_exact*ones(size(NN)), "--", NN, m_flow_vect,'-o', "LineWidth",3)
title("Comparison of numerical mass flow rate and exact value vs number of grid nodes \(N\)", "interpreter", "latex", FontSize=25)
xlabel("Number of nodes \(N\)", "Interpreter","latex", FontSize=20)
ylabel("Mass flow rate [kg/s]", "interpreter", "latex", FontSize=20)
grid on
legend("Exact value", "Numerical value",  "interpreter", "latex", fontsize=20);
set(gca,'FontSize',15)




