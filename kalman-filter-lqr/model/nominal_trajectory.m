%% 1.1 Trajectory with Jayachandran et al. model - theta_p estimate
% %% 1.1.1 Intialization
% clear
% bsa = 1.71;
% theta = [4.2;3.8;39.4;15.11;1;0.08;0.3287;8.2e9;0.4386;0.1207;0.5346;0.0782;84];
% % initial values (from Chapman et al.)
% % x0 = [0;0;0;1.2e11;1.2e11;1.2e11;1.2e11;theta(8)*(theta(7)/theta(10)-1)^(1/theta(9))];
% x0 = [0;0;0;1.2625e11;1.2625e11;1.2625e11;1.2625e11;0.2851e11];
% % x0 = [0;0;0;0;0;0;0;0];
% d = 50*bsa; % 6-MP dosage (mg 6-MP/m^2 body surface area)
% 
% %% 1.1.2 Calculation of trajectory for 5 cycles
% state = [];
% t = [];
% x0_i = x0;
% num_cycles = 34;
% num_t = 21; % number of timepoints to evaluate ODE45
% tspan = linspace(1,21,num_t);
% for i=1:num_cycles
%     [t_i,state_i] = ode45(@(t,x)jayachandran(t,x,d,theta),tspan,x0_i);
%     t_i = t_i + 20*(i-1); % shift the time index to the current cycle
%     % check lower bound
%     if state_i(end,8) < 1e9
%         d = 0.8*d;
%     elseif state_i(end,8)*0.6 > 2e9 % check upper bound
%         d = 1.2*d;
%     end
%     state = [state state_i];
%     t = [t t_i];
%     x0_i = state_i(end,:);
% end
% 
% clear i;
% 
% %% 1.1.3 Plot of trajectory for 5 cycles
% figure()
% hold on
% xlabel('Days Elapsed')
% ylabel('WBC Count')
% title("WBC Count vs Days Elapsed")
% for i=1:num_cycles
%     plot(t(:,i),state(:,8*i))
% end
% hold off
% 
% clear i;
% 
% figure()
% hold on
% xlabel('Days Elapsed')
% ylabel('Neutrophil Count')
% title("Neutrophil Lower Bound vs Days Elapsed")
% for i=1:num_cycles
%     plot(t(:,i),0.4*state(:,8*i))
% end
% hold off
% 
% clear i;

% %% Trajectory with theta_p estimate
% 
% % WBC trajectory over 5 cycles; lower bound is 
% [time1,state1] = ode45(@(t,x)jayachandran(t,x,d),[0 20],x0);
% 
% % check lower bound
% if state1(end,8) < 1e9
%     d = 0.8*d;
%     sprintf("less")
% elseif state1(end,8)*0.6 > 2e9 % check upper bound
%     d = 1.2*d;
% end
% [time2,state2] = ode45(@(t,x)jayachandran(t,x,d),[0 20],state1(end,:));
% time2 = time2+20;
% 
% % check lower bound
% if state2(end,8)*0.4 < 1e9
%     d = 0.8*d;
%     sprintf("less")
% elseif state2(end,8)*0.6 > 2e9 % check upper bound
%     d = 1.2*d;
% end
% [time3,state3] = ode45(@(t,x)jayachandran(t,x,d),[0 20],state2(end,:));
% time3 = time3+20*2;
% 
% % check lower bound
% if state3(end,8)*0.4 < 1e9
%     d = 0.8*d;
%     sprintf("less")
% elseif state3(end,8)*0.6 > 2e9 % check upper bound
%     d = 1.2*d;
% end
% [time4,state4] = ode45(@(t,x)jayachandran(t,x,d),[0 20],state3(end,:));
% time4 = time4+20*3;
% 
% % check lower bound
% if state4(end,8)*0.4 < 1e9
%     d = 0.8*d;
%     fprintf("less")
% elseif state4(end,8)*0.6 > 2e9 % check upper bound
%     d = 1.2*d;
% end
% [time5,state5] = ode45(@(t,x)jayachandran(t,x,d),[0 20],state4(end,:));
% time5 = time5+20*4;
% 
% % [time1_ext,state1_ext] = ode45(@(t,x)jayachandran(t,x),[0 200],x0);
% 
% figure()
% hold on
% plot(time1,state1(:,8))
% [state1_max,state1_max_idx] = max(state1);
% state1_max
% plot(time2,state2(:,8))
% plot(time3,state3(:,8))
% plot(time4,state4(:,8))
% plot(time5,state5(:,8))
% % plot(time1_ext,state1_ext(:,8))
% xlabel('Days Elapsed')
% ylabel('WBC Count')
% legend("cycle 1","cycle 2","cycle 3","cycle 4","cycle 5")
% title("WBC Count vs Days Elapsed")
% hold off
% 
% figure()
% hold on
% plot(time1,0.4*state1(:,8))
% plot(time2,0.4*state2(:,8))
% plot(time3,0.4*state3(:,8))
% plot(time4,0.4*state4(:,8))
% plot(time5,0.4*state5(:,8))
% xlabel('Days Elapsed')
% ylabel('ANC')
% legend("cycle 1","cycle 2","cycle 3","cycle 4","cycle 5")
% title("Neutrophil Lower Bound vs Days Elapsed")
% hold off
% 
% %% Moar Plots
% [time1,state1] = ode45(@(t,x)jayachandran(t,x,d),[0 200],x0);
% figure()
% hold on
% plot(time1,state1(:,1))
% plot(time1,state1(:,2))
% hold off
% 
% figure()
% hold on
% plot(time1,state1(:,3))
% hold off
% 
% figure()
% hold on
% plot(time1,state1(:,4))
% plot(time1,state1(:,5))
% plot(time1,state1(:,6))
% plot(time1,state1(:,7))
% hold off
% 
% figure()
% hold on
% plot(time1,state1(:,8))
% hold off
%% 1.2 Trajectory with Jayachandran et al. model - modified treatment protocol
%% 1.2.1 Intialization
clear
clc
bsa = 1.71;
theta = [4.2;3.8;39.4;15.11;1;0.08;0.3287;8.2e9;0.4386;0.1207;0.5346;0.0782;84];
% initial values (from Chapman et al.)
% x0 = [0;0;0;1.2e11;1.2e11;1.2e11;1.2e11;theta(8)*(theta(7)/theta(10)-1)^(1/theta(9))];
x0 = [0;0;0;1.2625e11;1.2625e11;1.2625e11;1.2625e11;0.2851e11];
% x0 = [0;0;0;0;0;0;0;0];
dose = 500000*bsa; % 6-MP dosage (mg 6-MP/m^2 body surface area)
%% 
%u = @(t,d) d*(t<=14) + 0*(t>14);
H = @(x) (x >= 0);
u = @(t,d) H(t)-H(t-14);
% define matrices and nonlinear components
A = [-theta(1)         0 0 0 0 0 0 0; ...
    theta(1) -theta(2)   0 0 0 0 0 0; ...
    0 0 -theta(6)          0 0 0 0 0; ...
    0 0 0 -theta(10)         0 0 0 0; ...
    0 0 0 theta(10) -theta(10) 0 0 0; ...
    0 0 0 0 theta(10) -theta(10) 0 0; ...
    0 0 0 0 0 theta(10) -theta(10) 0; ...
    0 0 0 0 0 0 theta(10) -theta(11)];

B = [5875.4;0;0;0;0;0;0;0];

f_hat = @(t,x) [0;...
    -theta(3)*x(2)/(theta(4)+x(2));...
    theta(5)*theta(3)*x(2)/(theta(4)+x(2));...
    (theta(7)*(theta(8)^theta(9))/(theta(8)^theta(9)+x(8)^theta(9))-theta(12)*x(3)/(theta(13)+x(3)))*x(4);...
    0;...
    0;...
    0;...
    0];

f = @(t,x,d) A*x + B*0.2*bsa*u(t) + f_hat(t,x);
%% 1.2.2 Calculation of trajectory for 5 cycles
state = [];
t = [];
x0_i = x0;
num_cycles = 7;
num_t = 21; % number of timepoints to evaluate ODE45
tspan = linspace(1,21,num_t);
for i=1:num_cycles
    [t_i,state_i] = ode45(@(t,x)f(t,x,dose),tspan,x0_i);
    t_i = t_i + 20*(i-1); % shift the time index to the current cycle
    % check lower bound
    if state_i(end,8) < 1.5e9
        dose = 0.8*dose;
    elseif state_i(end,8) > 3.0e9 % check upper bound
        dose = 1.2*dose;
    end
    x_dot = @(t,x)f(t,x,dose);
    state = [state state_i];
    t = [t t_i];
    x0_i = state_i(end,:);
end % sanity check code

clear i;

%% 1.2.3 Plot of trajectory for 5 cycles
figure()
hold on
xlabel('Days Elapsed')
ylabel('WBC Count')
title("WBC Count vs Days Elapsed")
for i=1:num_cycles
    plot(t(:,i),state(:,8*i))
end
hold off

clear i;

% figure()
% hold on
% xlabel('Days Elapsed')
% ylabel('Neutrophil Count')
% title("Neutrophil Lower Bound vs Days Elapsed")
% for i=1:num_cycles
%     plot(t(:,i),0.4*state(:,8*i))
% end
% hold off

clear i;

%% 2. Trajectory with Jost et al. model
%% 2.1 Intialization of constants
clear
bsa = 1.71;
theta = [31.2;12.72;0.019;9.9216;0.219*(bsa^1.16);2.06;0.146;0.103;0.866;2.3765];
% initial values (from Jost et al.)
x0 = [0;...
    0;...
    0;...
    theta(6)*theta(10)/theta(7);...
    theta(6)*theta(10)/theta(7);...
    theta(6)*theta(10)/theta(7);...
    theta(6)*theta(10)/theta(7);...
    theta(6)];

num_cycles = 34;

%% 2.2 Calculation of reference trajectory for 5 cycles

x_ref = [];
t_ref = [];
u_ref = [];
step_size_ref = 0.01; % very fine step size
num_t_ref = (1/step_size_ref)*21+1; % number of timepoints to evaluate ODE45
tspan = linspace(0,21,num_t_ref);

u_i = 50*bsa; % 6-MP dosage (mg 6-MP/m^2 body surface area)
x0_i = x0;

for i=1:num_cycles
    u_ref = [u_ref [u_i*(ones((1/step_size_ref)*14,1));zeros((1/step_size_ref)*(20-14)+1,1)]];
    [t_i,x_i] = ode45(@(t,x)jost(t,x,u_i,theta),tspan,x0_i);
    t_i = t_i + 21*(i-1); % shift the time index to the current cycle
    % check lower bound
    if x_i(end,8) < 1
        u_i = 0.8*u_i;
    elseif x_i(end,8) > 2 % check upper bound
        u_i = 1.2*u_i;
    end
    x_ref = [x_ref x_i];
    t_ref = [t_ref t_i];
    x0_i = x_i(end,:);
end

clear i t_i x_i x0_i u_i tspan;

% flatten to single nx8 matrix
x_ref_flattened = [];
t_ref_flattened = [];
u_ref_flattened = [];
for i=1:num_cycles
    x_ref_flattened = [x_ref_flattened; x_ref(:,8*(i-1)+1:8*i)];
    t_ref_flattened = [t_ref_flattened; t_ref(:,i)];
    u_ref_flattened = [u_ref_flattened; u_ref(:,i)];
end
clear i;

%% 2.3 Plot of 6-MP trajectory

figure()
hold on
xlabel('Days Elapsed')
ylabel('Amount (mg) or Concentration (mg/L blood)')
title("6-MP and 6-TGN Amount vs Days Elapsed")
plot(t_ref_flattened,x_ref_flattened(:,1))
plot(t_ref_flattened,x_ref_flattened(:,2))
plot(t_ref_flattened,x_ref_flattened(:,3))
legend("x_1 6-MP in the gut (mg)","x_2 6-MP in the bloodstream (mg)","x_3 6-TGN in the bloodstream (mg/L blood)")
hold off

%% 2.4 Plot of number of proliferating cells, number of cells in each compartment, and number of mature neutrophils

figure()
hold on
xlabel('Days Elapsed')
ylabel('Number of cells / L blood')
title("Number of Cells per L Blood vs Days Elapsed")
plot(t_ref_flattened,x_ref_flattened(:,4)*1e9)
plot(t_ref_flattened,x_ref_flattened(:,5)*1e9)
plot(t_ref_flattened,x_ref_flattened(:,6)*1e9)
plot(t_ref_flattened,x_ref_flattened(:,7)*1e9)
plot(t_ref_flattened,x_ref_flattened(:,8)*1e9)
legend("x_4 Proliferating cells","x_5 Compartment 1","x_6 Compartment 2","x_7 Compartment 3","x_8 Mature neutrophils")
hold off


%% 2.5 Plot of neutrophil trajectory (color-coded by cycle)

figure()
hold on
xlabel('Days Elapsed')
ylabel('ANC (Number of cells / L blood)')
title("All Neutrophil Count (ANC) x_8 vs Days Elapsed")
for i=1:num_cycles
    plot(t_ref(:,i),x_ref(:,8*i)*10^9)
end
yline(1e9,'-b','Desired lower bound')
yline(2e9,'-r','Desired upper bound')
hold off

clear i;

%% 2.6 Linearization

bsa = 1.71; % trajectory near original bsa
theta = [31.2;12.72;0.019;9.9216;0.219*(bsa^1.16);2.06;0.146;0.103;0.866;2.3765];
x_lin = [];
t_lin = [];
u_lin = [];
step_size_lin = 0.01; % coarse step size; NOTE: because of forward euler's small stability region, the largest I could get this is 0.05 with a different bsa
num_t_lin = (1/step_size_lin)*21+1; % number of timepoints to evaluate ODE45
tspan = linspace(0,21,num_t_lin);


t_lin_j_start = 0; % current time at end of last cycle
x_t = x0; % initial "guess"
u_i = 50*bsa;


dx_t = zeros(8,1); % keep track of delta x
%du_t = 0; % keep track of delta u

dx_t_arr = [];

for i=1:num_cycles
    u_lin = [u_lin u_i];
    u_star = u_ref(end,i);
    du_t = u_i - u_star;

    x_ref_i = x_ref(:,8*(i-1)+1:8*i);
    x_lin_j = [];
    t_lin_j = [];
    
    for j=1:length(tspan)
        x_star = transpose(x_ref_i(round((j-1)*step_size_lin/step_size_ref)+1,:));
        %round((j-1)*step_size_lin/step_size_ref)+1
        x_t = x_star+dx_t;
        x_lin_j = [x_lin_j; transpose(x_t)];

        [f, dfdx, dfdu, dgdx] = jost_fwd_euler(tspan(j),x_star,u_star,theta,step_size_lin);
        dx_t = dfdx*dx_t+dfdu*du_t;
        if tspan(j) >= 14
            % u_i and u_ref(i) are both zero after 14 days
            du_t = 0;
        end
        t_lin_j = [t_lin_j; t_lin_j_start+tspan(j)];
    end
    t_lin = [t_lin t_lin_j];
    x_lin = [x_lin x_lin_j];

    % check lower bound
    if x_t(8) < 1
        u_i = 0.8*u_i;
    elseif x_t(8) > 2 % check upper bound
        u_i = 1.2*u_i;
    end
    t_lin_j_start = t_lin_j_start + 21; % shift the time index to the current cycle
end

t_lin_flattened = [];
x_lin_flattened = [];
for i=1:num_cycles
    x_lin_flattened = [x_lin_flattened; x_lin(:,8*(i-1)+1:8*i)];
    t_lin_flattened = [t_lin_flattened; t_lin(:,i)];
end

figure()
hold on
xlabel('Days Elapsed')
ylabel('Number of cells / L blood')
title("Number of Cells per L Blood vs Days Elapsed")
plot(t_lin_flattened,x_lin_flattened(:,4)*1e9)
plot(t_lin_flattened,x_lin_flattened(:,5)*1e9)
plot(t_lin_flattened,x_lin_flattened(:,6)*1e9)
plot(t_lin_flattened,x_lin_flattened(:,7)*1e9)
plot(t_lin_flattened,x_lin_flattened(:,8)*1e9)
legend("x_4 Proliferating cells","x_5 Compartment 1","x_6 Compartment 2","x_7 Compartment 3","x_8 Mature neutrophils")
hold off


figure()
hold on
xlabel('Days Elapsed')
ylabel('ANC (Number of cells / L blood)')
title("All Neutrophil Count (ANC) x_8 vs Days Elapsed")
plot(t_lin_flattened,x_lin_flattened(:,8)*1e9)
plot(t_ref_flattened,x_ref_flattened(:,8)*1e9)

yline(1e9,'-b')
yline(2e9,'-r')

legend("Linearized","Original","Desired lower bound","Desired upper bound","Location","southwest")
hold off

clear i j x_lin_j x0_i u_i tspan;

%% 2.7.1 Trajectory with Noise; Setup

num_cycles = 14; %34

step_size_noisy = step_size_ref;
num_t_noisy = (1/step_size_noisy)*21+1; % number of timepoints to evaluate ODE45
tspan = linspace(0,21,num_t_noisy);

var_w = 0.01; % TODO: change this
var_v = 0.774341747; % CHECK if documented

%% 2.7.2 Trajectory with Noise; Reactive Controller
tic

u_i = 50*bsa;
x0_i = x0;

y = [];

x_noisy_r = [];
t_noisy_r = [];
u_noisy_r = [];


for i=1:num_cycles
    u_i_all = [transpose(repelem(u_i,(1/step_size_noisy)*14));transpose(repelem(0,(1/step_size_noisy)*(21-14)+1))];
    u_noisy_r = [u_noisy_r u_i_all];
    [t_i,x_i] = ode45(@(t,x)jost_noisy(t,x,u_i,theta,step_size_noisy,var_w),tspan,x0_i);
    t_i = t_i + 20*(i-1); % shift the time index to the current cycle
    
    % make a "reading"
    % ground truth (at the end of the cycle)
    y_true = x_ref(end,8*i);
    % sample a reading
    y_i = normrnd(y_true,var_v);
    y = [y y_i];

    % check lower bound
    if y_i < 1
        u_i = 0.8*u_i;
    elseif y_i > 2 % check upper bound
        u_i = 1.2*u_i;
    end
    x_noisy_r = [x_noisy_r x_i];
    t_noisy_r = [t_noisy_r t_i];
    
    x0_i = x_i(end,:);
end

clear i t_i x_i x0_i u_i tspan;

% flatten to single nx8 matrix
x_noisy_r_flattened = [];
t_noisy_r_flattened = [];
u_noisy_r_flattened = [];
for i=1:num_cycles
    x_noisy_r_flattened = [x_noisy_r_flattened; x_noisy_r(:,8*(i-1)+1:8*i)];
    t_noisy_r_flattened = [t_noisy_r_flattened; t_noisy_r(:,i)];
    u_noisy_r_flattened = [u_noisy_r_flattened; u_noisy_r(:,i)];
end
clear i;

toc

%%
figure()
hold on
xlabel('Days Elapsed')
ylabel('ANC (Number of cells / L blood)')
title("All Neutrophil Count (ANC) x_8 vs Days Elapsed")
plot(t_noisy_r_flattened,x_noisy_r_flattened(:,8)*1e9)
plot(t_ref_flattened,x_ref_flattened(:,8)*1e9)
yline(1e9,'-b')
yline(2e9,'-r')
legend("Noisy","Original","Desired lower bound","Desired upper bound","Location","southwest")
hold off

figure()
hold on
title("Noisy model reactive controller control sequence")
plot(t_noisy_r_flattened,u_noisy_r_flattened)
xlabel("Time elapsed (days)")
ylabel("Input dosage (mg)")
hold off

%% 2.7.3 Trajectory with Noise; KF Controller
tic

% Simu_symlation Length
t_s = step_size_noisy;
t_end = 21*num_cycles;
N = t_end/t_s;
time = [0:t_s:t_end];

% Dimension
% n: dimension of state vector x
n = 8;
% m: dimension of observation vector y
m = 1;

% Initialize motion function
syms x1_sym x2_sym x3_sym x4_sym x5_sym x6_sym x7_sym x8_sym u_sym

f = [-theta(1)                0 0 0 0 0 0 0; ...
    theta(1) -theta(2)          0 0 0 0 0 0; ...
    0 theta(3)*theta(4) -theta(5) 0 0 0 0 0; ...
    0 0 0 -theta(7)                 0 0 0 0; ...
    0 0 0 theta(7) -theta(7)          0 0 0; ...
    0 0 0 0 theta(7) -theta(7)          0 0; ...
    0 0 0 0 0 theta(7) -theta(7)          0; ...
    0 0 0 0 0 0 theta(7) -theta(10)]*[x1_sym;x2_sym;x3_sym;x4_sym;x5_sym;x6_sym;x7_sym;x8_sym] + ...
    [0.22;0;0;0;0;0;0;0]*u_sym + ...
    [0;...
    0;...
    0;...
    theta(7)*((theta(6)/x8_sym)^theta(9))*x4_sym-theta(7)*theta(8)*((theta(6)/x8_sym)^theta(9))*x3_sym*x4_sym;...
    0;...
    0;...
    0;...
    0];

f = [x1_sym;x2_sym;x3_sym;x4_sym;x5_sym;x6_sym;x7_sym;x8_sym] + t_s * f;

dfdx = jacobian(f,[x1_sym;x2_sym;x3_sym;x4_sym;x5_sym;x6_sym;x7_sym;x8_sym]);
dfdu = jacobian(f,u_sym);

% Initialize observation function
g = [x8_sym];
dgdx = jacobian(g,[x1_sym;x2_sym;x3_sym;x4_sym;x5_sym;x6_sym;x7_sym;x8_sym]);

% Initialize disturbance xi and covariance Q
xi = zeros([n 1 N]);
Q = zeros([n n N]);
Q_const = var_w * eye(n);
for iter = 1:N
    Q(:,:,iter) = Q_const;
    xi(:,:,iter) = Q(:,:,iter) * randn([n,1]);
end
clear Q_const

% Initialize measurement noise n and covariance R
nn = zeros([m 1 N]);
R = zeros([m m N]);
R_const = var_v * eye(m);
for iter = 1:N
    R(:,:,iter) = R_const;
    nn(:,:,iter) = R(:,:,iter) * randn([m,1]);
end
clear R_const

% Initialize nominal trajectory and control to track
% nominal control is set as linear for test purposes
A_nominal = zeros([n n N]);
b_nominal = zeros([n N]);

% convert and reshape for Sam's code
u_nominal = transpose(u_ref_flattened(1:N,1));

x_nominal = zeros([n 1 N]);
% Simulate System for the nominals
x_o = x0;
for k = 1:N
    % Compute linearized c2d-ed state matrices, we use this to simulate the
    % system because we do not know the discrete time nonlinear equation
    % for the inverse pendulum problem
    if k == 1
        x_tmp = x_o;
    else
        x_tmp = x_nominal(:,:,k-1);
    end
    [A_actual,b_actual,~] = c2d_jost(x1_sym,x2_sym,x3_sym,x4_sym,x5_sym,x6_sym,x7_sym,x8_sym,u_sym,dfdx,dfdu,dgdx,x_tmp,u_nominal(:,k),t_s);
    A_nominal(:,:,k) = A_actual;
    b_nominal(:,k) = b_actual;
    % Time step
    x_nominal(:,:,k) = double(subs(f,{x1_sym,x2_sym,x3_sym,x4_sym,x5_sym,x6_sym,x7_sym,x8_sym,u_sym},{x_tmp(1),x_tmp(2),x_tmp(3),x_tmp(4),x_tmp(5),x_tmp(6),x_tmp(7),x_tmp(8),u_nominal(:,k)}));
end
clear k x_tmp A_actual b_actual


% Initialize Cost Parameters
W = zeros([n n N+1]);
rho = zeros([n 1 N+1]);
lambda = zeros([1 N]);
W_const = 100 * eye(n);
lambda_const = 0.1;
for iter = 1:N+1
    W(:,:,iter) = W_const;
    if iter == 1
        rho(:,:,iter) = x_o;
    else
        rho(:,:,iter) = x_nominal(:,:,iter-1);
    end
end
for iter = 1:N
    lambda(:,iter) = lambda_const;
end
clear W_const lambda_const

% Note: Iteration variables does not correspond to the notes due to matlab
% indexing start with 1 instead of 0

x_o = x0;
x_hat_o = x0; 
sigma_o = eye(n);

% Simulate System with control
is_controlled = 1;
[x_c,y_c,x_hat_p_c,u_c] = simulate_jost(N,n,m,t_s,x1_sym,x2_sym,x3_sym,x4_sym,x5_sym,x6_sym,x7_sym,x8_sym,u_sym,f,dfdx,dfdu,dgdx,xi,Q,nn,R,W,rho,lambda,x_o,x_hat_o,sigma_o,A_nominal,b_nominal,x_nominal,u_nominal,is_controlled);


% Process Result
x1_c_plot = [x_o(1,1);reshape(x_c(1,1,:),[N 1])];
x2_c_plot = [x_o(2,1);reshape(x_c(2,1,:),[N 1])];
x3_c_plot = [x_o(3,1);reshape(x_c(3,1,:),[N 1])];
x4_c_plot = [x_o(4,1);reshape(x_c(4,1,:),[N 1])];
x5_c_plot = [x_o(5,1);reshape(x_c(5,1,:),[N 1])];
x6_c_plot = [x_o(6,1);reshape(x_c(6,1,:),[N 1])];
x7_c_plot = [x_o(7,1);reshape(x_c(7,1,:),[N 1])];
x8_c_plot = [x_o(8,1);reshape(x_c(8,1,:),[N 1])];
y1_c_plot = [x_o(1,1);reshape(y_c(1,1,:),[N 1])];
x1_hat_c_plot = [x_hat_o(1,1);reshape(x_hat_p_c(1,1,:),[N 1])];
x2_hat_c_plot = [x_hat_o(2,1);reshape(x_hat_p_c(2,1,:),[N 1])];
x3_hat_c_plot = [x_hat_o(3,1);reshape(x_hat_p_c(3,1,:),[N 1])];
x4_hat_c_plot = [x_hat_o(4,1);reshape(x_hat_p_c(4,1,:),[N 1])];
x5_hat_c_plot = [x_hat_o(5,1);reshape(x_hat_p_c(5,1,:),[N 1])];
x6_hat_c_plot = [x_hat_o(6,1);reshape(x_hat_p_c(6,1,:),[N 1])];
x7_hat_c_plot = [x_hat_o(7,1);reshape(x_hat_p_c(7,1,:),[N 1])];
x8_hat_c_plot = [x_hat_o(8,1);reshape(x_hat_p_c(8,1,:),[N 1])];

toc
%% Plotting 2.7.3

figure()
plot(t_noisy_r_flattened(1:N),x_noisy_r_flattened(1:N,8));
hold on
rho8_plot = reshape(rho(8,1,:),[N+1 1]);
plot(time,x8_c_plot)
hold on
plot(time,rho8_plot, "LineWidth", 3, "Color", "#7E2F8E")
xlabel('Days Elapsed')
ylabel('ANC (Number of cells / L blood)')
title("All Neutrophil Count (ANC) Comparison")
legend('Reactive','Kalman Filter','Nominal')

figure()
hold on
title("Noisy Model KF Controller Control Sequence")
plot(time(1:N),u_c)
xlabel("Time elapsed (days)")
ylabel("Input dosage (mg)")
hold off

figure()
hold on
title("Noisy Model Control Sequence Comparison")
plot(t_noisy_r_flattened,u_noisy_r_flattened)
plot(time(1:N),u_c)
legend('reactive control sequence','KF control sequence')
xlabel("Time elapsed (days)")
ylabel("Input dosage (mg)")
hold off

%% 2.8.1 Control Input Smoothing - Cyclic Intervals

N = length(u_c);
u_smoothed_flattened = zeros(1, N);
u_smoothed = [];
t = time(1:N);

i=1;
start_i = 1;
day = 14;
t_i = t(i);
sum = 0;
n_points = 0;

while i <= N
    sum = sum + u_c(i);
    if t(i) == day
        u_i = sum/n_points;
        u_smoothed_flattened(start_i:i) = u_i;
        idx_7_days_later = find(t == day+6);
        u_smoothed_flattened(i:idx_7_days_later) = 0;
        u_smoothed = [u_smoothed u_i];

        i = idx_7_days_later;
        n_points = 0;
        sum = 0;
        start_i = i+1;
        day = day + 20;
        if day > N
            day = N;
        end
    else
        n_points = n_points + 1;
    end
    i = i + 1;
end

if n_points > 0
    u_smoothed_flattened(start_i:i-1) = sum/n_points;
    u_smoothed = [u_smoothed u_i];
end

figure()
hold on
plot(time(1:N),u_c);
plot(t,u_smoothed_flattened);
lgd = legend('KF control sequence', 'Smoothed KF control sequence');
lgd.Location = "northoutside";
xlabel("Time elapsed (days)");
ylabel("Input dosage (mg)");
title("Smoothed Control Sequence for Kalman Filter Controller")
hold off

figure()
hold on
plot(t_noisy_r_flattened,u_noisy_r_flattened, "LineWidth", 2)
plot(t,u_smoothed_flattened, "LineWidth", 2)
lgd = legend('Reactive control sequence', 'Smoothed KF control sequence');
lgd.Location = "northoutside";
xlabel("Time elapsed (days)");
ylabel("Input dosage (mg)");
title("Smoothed KF Control Sequence vs. Reactive Control Sequence")
hold off

clear start_i i day idx_7_days_later sum u_i


%% 2.8.2 Model's Response with Smoothed Inputs - Cyclic Intervals

x_smoothed = [];
t_smoothed = [];
step_size_smoothed = 0.01; % very fine step size
num_t_smoothed = (1/step_size_smoothed)*21+1; % number of timepoints to evaluate ODE45
tspan = linspace(0,21,num_t_smoothed);

x0_i = x0;
    
for i=1:num_cycles
    [t_i,x_i] = ode45(@(t,x)jost(t,x,u_smoothed(i),theta),tspan,x0_i);
    t_i = t_i + 21*(i-1); % shift the time index to the current cycle

    x_smoothed = [x_smoothed x_i];
    t_smoothed = [t_smoothed t_i];
    x0_i = x_i(end,:);
end

clear i t_i x_i x0_i u_i tspan;

% flatten to single nx8 matrix
x_smoothed_flattened = [];
t_smoothed_flattened = [];
for i=1:num_cycles
    x_smoothed_flattened = [x_smoothed_flattened; x_smoothed(:,8*(i-1)+1:8*i)];
    t_smoothed_flattened = [t_smoothed_flattened; t_smoothed(:,i)];
end
clear i;

%% 2.8.3 Plotting Response of Smoothed Inputs

figure()
hold on
xlabel('Days Elapsed')
ylabel('ANC (Number of cells / L blood)')
title("All Neutrophil Count (ANC) x_8 vs Days Elapsed")
% plot(time,rho8_plot*1e9, "LineWidth", 2)
plot(time,x8_c_plot*1e9);
% plot(t_noisy_r_flattened,x_noisy_r_flattened(:,8)*1e9)
plot(t_smoothed_flattened,x_smoothed_flattened(:,8)*1e9, "LineWidth", 2)
yline(1e9,'-b')
yline(2e9,'-r')
legend("Original", "Smoothed", "Desired lower bound","Desired upper bound","Location","southwest")
hold off

figure()
hold on
xlabel('Days Elapsed')
ylabel('ANC (Number of cells / L blood)')
title("All Neutrophil Count (ANC) x_8 vs Days Elapsed")
plot(time,rho8_plot*1e9, "LineWidth", 2)
% plot(t_noisy_r_flattened,x_noisy_r_flattened(:,8)*1e9)
plot(t_smoothed_flattened,x_smoothed_flattened(:,8)*1e9, "LineWidth", 2)
yline(1e9,'-b')
yline(2e9,'-r')
legend("Nominal", "Smoothed", "Desired lower bound","Desired upper bound","Location","southwest")
hold off



