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
alpha_opt = fminsearch(@nozzle1d_alpha, [0.1,0.13]);  %alpha_u, alpha_p

fprintf("The optimal value according to Nelder Mead (check the initial guess) alpha_p is %.5f, while the optimal value for alpha_u is %.5f\n", alpha_opt(2), alpha_opt(1));


 
alpha_p_vect = linspace(0.1,0.25, 10);    %it s the x-axis
alpha_u_vect = linspace(0.1, 1.7, 10);



tic   %69 secondi circa sul mio pc per un 100x100
for i =1:length(alpha_u_vect)
    for j =1:length(alpha_p_vect)
    
   
    alpha_opt = [alpha_u_vect(i), alpha_p_vect(j)];
    it = nozzle1d_alpha(alpha_opt);
    it_matrix(i,j) = it;
    end
end


toc
[U,P] = meshgrid(alpha_u_vect, alpha_p_vect);
figure
surf(U,P,it_matrix');
colormap winter
ylabel("\(\alpha_p\)", "Interpreter","latex", FontSize=20)
xlabel("\(\alpha_u\)", "interpreter", "latex", fontsize=20)
zlabel("Number of iterations to reach a tolerance of \(1e-6\)", "Interpreter","latex", FontSize=20)
title("Required iterations vs under-relaxation coefficients, \(N=21\), \(tol = 1e-6\)", "Interpreter","latex", FontSize=20)
set(gca,'FontSize',15)
hold on
min_iter = min(min(it_matrix));
[pos_u,pos_p]=find(it_matrix == min_iter);

alpha_u_opt = alpha_u_vect(pos_u);
alpha_p_opt = alpha_p_vect(pos_p);

fprintf("The grid search yields alpha_u_opt = %.5f, alpha_p_opt = %.5f\n", alpha_u_opt(1), alpha_p_opt(1))
fprintf("The miminum found is %d iterations\n", min(min(it_matrix)))

plot3(alpha_u_opt(1),alpha_p_opt(1), min(min(it_matrix)), '.r','markersize',50)


P = it_matrix;

for i=1:size(P,1)
    for j = 1:size(P,2)
        if (P(i,j) >1000)
            P(i,j) = P(i,j)-1000;
        end
    end
end

disp(sum(sum(P)))





