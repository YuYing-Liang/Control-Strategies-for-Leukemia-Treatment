%% 1.1 Adaptive Control with Least Squares (linear case)

% Initialization
clear all;
clc; 
close all; 
format longE

N = 15;
time = (0:1:N);

n = 3; % dimension of state vector x 
m = 3; % dimension of observation vector y

p = [-5;15];

% 1. Derive Matrices for system (A,B and C) - constant for N timesteps
dt = 1;
A_c = [
    p(1)  dt  p(2)*dt^2; 
    0     p(1)  dt; 
    0     0     p(1);
];
A_actual = [
    1   dt  0.5*dt^2; 
    0   1   dt; 
    0   0   1;
];
b_c = [0;0;1];
C_c = eye(m,n);

% 2. Generate disturbance and measurement noise given covariances Q,R

sigQ = 1;
sigR = 5;
Q_const = sigQ * eye(n);
R_const = sigR * eye(m);

d = zeros([n 1 N]);                 % disturbance
mn = zeros([m 1 N]);                % measurement noise

Q = repmat(Q_const, [1 1 N+1]);
R = repmat(R_const, [1 1 N+1]);

for i = 1:N
    d(:,:,i) = Q_const * randn([n,1]);
    mn(:,:,i) = R_const * randn([m,1]);
end
clear sigQ sigR Q_const R_const i

% 3. Generate target trajectory
rho = zeros([n 1 N+1]);

a_t = -15;                          % target acceleration
sp = [100;0;-10];                   % starting point [x, x_dot, x_dot_dot]
for i = 1:N+1
    t = (i - 1) * dt;
    rho(:,:,i) = [
        sp(1,1) + sp(2,1)*t + 0.5*(t^2)*a_t;
        sp(2,1) + t*a_t;
        a_t;
    ];
end
clear i

% 4. Simulate System
x_calc = zeros([n,1,N+1]);
y_calc = zeros([m,1,N+1]);
p_calc = zeros([length(p),1,N+1]);
u_calc = zeros([1,N]);

x_calc(:,:,1) = sp;
y_calc(:,:,1) = sp;
p_calc(:,:,1) = p;

for k=1:N     
    % find vector p_star that minimizes measurement error
    A_guess = @(x)[
        x(1)   dt  x(2)*dt^2; 
        0      x(1)  dt; 
        0      0     x(1);
    ];
    f = @(x)calc_p_star(A_guess(x),A_actual,b_c,C_c,sp,mn,u_calc,k);
    [p_star, fval] = fmincon(f,p_calc(:,:,k));
    
    % find input u that minimizes output error
    A_star = [
        p_star(1)   dt          p_star(2)*dt^2; 
        0           p_star(1)   dt; 
        0           0           p_star(1);
    ];
    g = @(u)calc_adaptive_control(A_star, b_c, C_c, rho, x_calc(:,:,k), u, 10, 0.001, 5, 5, k);
    [u_star, gval] = fmincon(g,repmat(u_calc(k), [5 1]));
    
    u_calc(k+1) = u_star(1);
    p_calc(:,:,k+1) = p_star;
    
    % this is suppose to be the true system that the controller doesn't
    % know - controller is trying to learn the underlying true system
    x_calc(:,:,k+1) = A_actual*x_calc(:,:,k) + b_c*u_calc(k+1);
    y_calc(:,:,k+1) = C_c*x_calc(:,:,k+1);  
end

%% 1.2 Plots
x1 = reshape(x_calc(1,1,:),[1 N+1]);
x2 = reshape(x_calc(2,1,:),[1 N+1]);
x3 = reshape(x_calc(3,1,:),[1 N+1]);

figure(1)
rho1 = reshape(rho(1,1,:),[1 N+1]);
plot(time,x1)
hold on
plot(time,rho1)
title('Position over Time')
xlabel('Time t_k')
ylabel('Distance (m)')
legend('adaptive MPC','target')
hold off

p1 = reshape(p_calc(1,1,:),[1 N+1]);
p2 = reshape(p_calc(2,1,:),[1 N+1]);
figure(2)
plot(time,p1)
hold on
plot(time,p2)
title('Model Parameters')
xlabel('Time t_k')
ylabel('Paramter Value')
legend('Expected to be 1','Expected to be 0.5')
hold off

figure(3)
plot(time,x3)
hold on
plot(time,u_calc)
title('Acceleration over Time')
xlabel('Time t_k')
ylabel('acceleration (m/s^2)')
legend('x_3 - Actual acceleration','u - Control Input')
hold off

%% 2.1 Adaptive control for a Non-linear system
