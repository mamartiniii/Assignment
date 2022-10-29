close all
clear 
clc

%% Data

L=2;       %[m]
rho=1;     %[kg/m^3]
A_in=0.5;  %[m^2]
A_out=0.1; %[m^2]
p_0=10;    %[Pa]
p_out=0;   %[Pa]  NB: pressure is relative

%% Grid
N=21;  % 100 is a perfect value
n=N-1;
x_p=linspace(0,L,N);
x_u=linspace(x_p(2)/2,L-x_p(2)/2,n);
x_p=x_p';
x_u=x_u';


A_p=zeros(N,1);
A_u=zeros(n,1);
A=@(x) A_in + ((A_out- A_in)/L * x);
A_p=A(x_p);
A_u=A(x_u);

%% Initialize velocity and pressure
u_old=zeros(n,1);
p=zeros(N,1);

m=1;   %[kg/s]   initial guess for mass flow rate
u_old=m./(rho*A_u);
p=p_0-(p_0-p_out).*x_p/L;

%% solve the problem
alpha_p=0.1;   % under-relaxation coefficient-> to be tuned. 0.02 is a perfect value
alpha_u=0.1;   % under-relaxation coefficient-> to be tuned. 0.02 is a perfect value
it=0;
toll_u=1e-6;
toll_p=1e-6;
it_max=1000;
r_u=10;
r_p=10;
r_p_vect = [];
r_u_vect = [];

while ((r_u>toll_u && r_p>toll_p) && (it<it_max))  % controllare sta condizione 

    % Solve system Mu*u=b
    M_u=zeros(n);
    b_u=zeros(n,1);

   
    for i=2:n-1
        M_u(i,i)=rho*(u_old(i)+u_old(i+1))/2*A_p(i+1);
        M_u(i,i-1)=-rho*(u_old(i-1)+u_old(i))/2*A_p(i);
        b_u(i)=(p(i)-p(i+1))*A_u(i);
    end

    %CV1
    F_a=rho*u_old(1)*A_u(1);
    F_b=rho*(u_old(1)+u_old(2))/2*A_p(2);
    M_u(1,1)=F_b+F_a/2*(A_u(1)/A_p(1))^2;
    b_u(1)=(p_0-p(2))*A_u(1)+F_a*A_u(1)/A_p(1)*u_old(1);

    %CVn
    M_u(n,n)=rho*u_old(n)*A_u(n);
    M_u(n,n-1)=-rho*(u_old(n-1)+u_old(n))/2*A_p(N-1);
    b_u(n)=(p(N-1)-p(N))*A_u(n);

    u=M_u\b_u;


    % Solve system Mp*p'=b
    M_p=zeros(N);
    b_p=zeros(N,1);
    d_u=zeros(n,1);

    for i=1:n   %need a bit of explanation
        d_u(i)=A_u(i)/M_u(i,i);
    end

    for i=2:N-1
        a_W=rho*d_u(i-1)*A_u(i-1);
        a_E=rho*d_u(i)*A_u(i);
        M_p(i,i)=a_W+a_E;
        M_p(i,i-1)=-a_W;
        M_p(i,i+1)=-a_E;
        b_p(i)=rho*(A_u(i-1)*u(i-1)-A_u(i)*u(i));
    end

    %CV1
    M_p(1,1)=1;
    b_p(1)=0;

    %CVN
    M_p(N,N)=1;
    b_p(N)=0;

    p_first=M_p\b_p;

    % correct pressure and velocity
    u_calc=zeros(n,1);
    p_calc=zeros(N,1);
    for i=1:n
        u_calc(i)=u(i)+d_u(i)*(p_first(i)-p_first(i+1));
        if i>1
            p_calc(i)=p(i)+p_first(i);
        else
            p_calc(1)=p_0-rho/2*(u_calc(1)*A_u(1)/A_p(1))^2;
        end
    end

    %calculate relative residuals

    r_u=norm(M_u*u_calc-b_u)/norm(diag(M_u).*u_calc);
    r_p=norm(b_p);
    
    r_u_vect = [r_u_vect; r_u];
    r_p_vect = [r_p_vect; r_p];
    %update the values
    u_new=alpha_u*u_calc+(1-alpha_u)*u_old;
    p_new=alpha_p*p_calc+(1-alpha_p)*p;

    u_old=u_new;
    p=p_new;

    it=it+1;

end


% %% plot
% 
% figure
% plot(x_p,p_new,'Linewidth',2)
% title('pressure')
% 
% figure
% plot(x_u,u_new,'Linewidth',2)
% title('velocity')

%% Print relevant values

mean_m = mean(rho*u_new.*A_u);
m_in = rho*u_new(1)*A_in;
m_out = rho*u_new(end)*A_out;
check = (m_in - m_out)/mean_m;
fprintf("The final global mass balance is %.2f\n", check);
fprintf("Number of nodes: %d\n", N);
fprintf("Under-relaxation coefficients: \\alpha u = %.2f, \\alpha p = %.2f\n ", alpha_u, alpha_p)


%% Exact solution and plot

u_exact= A_out./A_u*sqrt(2*p_0/rho);
p_exact= p_0.*(1-(A_out./A_p).^2);

figure
subplot(1,2,1)
plot(x_u,u_new,x_u,u_exact,'Linewidth',3)
title('Velocity', "Interpreter","latex", FontSize=20)
legend('Numerical solution','Exact solution', fontsize=20)
xlabel('x [m]', "Interpreter","latex", FontSize=20)
ylabel('v [m/s]', "Interpreter","latex", FontSize=20)
grid on
set(gca,'FontSize',15)

subplot(1,2,2)
plot(x_p,p_new,x_p,p_exact,'Linewidth',3)
title('Pressure', "Interpreter","latex", Fontsize=20)
legend('Numerical solution','Exact solution', fontsize=20)
xlabel('x [m]', "Interpreter","latex", FontSize=20)
ylabel('P [Pa]', "Interpreter","latex", fontsize=20)
grid on
set(gca,'FontSize',15)



%% Plot residuals

figure;

subplot(1,2,1)

plot(r_u_vect, "Linewidth", 3);
title("Residual of momentum equation", "interpreter", "latex", fontsize=20)
grid on
xlabel("Iteration", "Interpreter","latex", fontsize=20);
ylabel("Residual \(r_u\)", "Interpreter","latex", FontSize=20);
set(gca,'FontSize',15)


subplot(1,2,2)
plot(r_p_vect, "LineWidth",3);
title("RHS of pressure correction equation", "interpreter", "latex", FontSize=20);
grid on
xlabel("Iteration", "Interpreter","latex", FontSize=20);
ylabel("RHS \(r_p\)", "Interpreter","latex", FontSize=20);
set(gca,'FontSize',15)


%% Relative errors for final u, p and and flow rate

% relative error for final u
err_u = abs(u_new(end) - u_exact(end));
err_u = err_u / u_exact(end);

% relative error for final p (too small numbers, it should be zero)
err_p = abs(p_new(end) - p_exact(end));
err_p = err_p / abs(p_exact(end));

% relative error for final flow rate
m_exact = A_out * sqrt(2 * rho * p_0);
err_m = abs(m_out - m_exact);
err_m = err_m / m_exact;


