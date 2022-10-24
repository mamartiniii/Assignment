function [u, p, x_u, x_p] = nozzle_1d(A, N, alpha_p, alpha_u, toll_u, toll_p, it_max, m)
%
% [u, p, x_u, x_p] = nozzle_1d(A, N, alpha_p, alpha_u, toll_u, toll_p, it_max, m)
%
% The function solves the problem of an incompressible flow in a converging
% nozzle. The SIMPLE algorithm is applied. This is a 1st order algorithm.
% INPUT
%
%
% OUTPUT
%
% By Mattia Martini, Francesco Sala, Nicolò Viscusi


% Given data (do not change)
L=2;       %[m]
rho=1;     %[kg/m^3]
p_0=10;    %[Pa]
p_out=0;   %[Pa]  NB: pressure is relative
n = N-1;
% Assembly of vector of nodes
x_p=linspace(0,L,N);
x_u=linspace(x_p(2)/2,L-x_p(2)/2,n);
x_p=x_p';
x_u=x_u';


%A_p=zeros(N,1);
%A_u=zeros(n,1);

A_p=A(x_p);
A_u=A(x_u);
%% Initialize velocity and pressure
%u_old=zeros(n,1);
%p=zeros(N,1);

   %[kg/s]   initial guess for mass flow rate
u_old=m./(rho*A_u);
p=p_0-(p_0-p_out).*x_p/L;

%% solve the problem

it=0;

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