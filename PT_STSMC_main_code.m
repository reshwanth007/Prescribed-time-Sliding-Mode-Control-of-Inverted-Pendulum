%% 1. WORKSPACE INITIALIZATION
clear; clc; close all;

%% 2. PARAMETER DEFINITIONS
% --- Plant Parameters ---
p.l = 1; p.k = 0.5; p.g = 9.81; p.m = 0.5;
I = p.m*p.l^2; % Moment of inertia

% --- PT-SMC Gains ---
c.c1 = 2.0;  c.c2 = 1.0;
c.c3 = 3.0;  c.c4 = 1.0;

% --- Super-Twisting Algorithm (STA) Gains ---
c.lambda = 10.0; % STA Proportional Gain
c.W = 5.0;       % STA Integral Gain

% --- Prescribed Times ---
T_as = 0.4;      % Prescribed time for sliding surface to converge (s)
T_a  = 0.5;      % Prescribed time for states (error) to converge (s)

% --- Simulation Parameters ---
T_sim = 10; dt = 1e-3; % Increased simulation time to 10s
N = round(T_sim/dt); t = (0:N-1)'*dt;

% --- Time-Varying Setpoint Definition ---
ref1_deg = 85;
ref2_deg = 170;

% --- Pre-allocation and Initialization ---
theta = zeros(N,1); theta_dot = zeros(N,1); tau = zeros(N,1);
setpoint_log_deg = zeros(N,1); % To log the setpoint for plotting
w = 0; % STA integral term initialization

%% 3. SIMULATION EXECUTION
for k = 2:N
    % --- Time-Varying Setpoint ---
    if t(k) < 5
        setpoint_deg = ref1_deg;
    else
        setpoint_deg = ref2_deg;
    end
    setpoint = setpoint_deg * pi/180;
    setpoint_log_deg(k) = setpoint_deg; % Log for plotting

    % --- State Errors ---
    e  = setpoint - theta(k-1);
    ed = 0 - theta_dot(k-1); % Reference velocity is zero

    % --- Time-Varying Gains Calculation (Using Original Functions) ---
    phi_a     = phi_fun(t(k), 0, T_a, 1);
    phi_a_dot = phi_dot_fun(t(k), 0, T_a, 1);
    phi_as    = phi_fun(t(k), 0, T_as, 1);

    % --- Sliding Surface ---
    S = ed + (c.c1 + c.c2*phi_a)*e;

    % --- Super-Twisting Algorithm ---
    u_sta = c.lambda * sqrt(abs(S)) * sign(S) + w;
    w_dot = c.W * sign(S);
    w = w + dt * w_dot; % Integrate w

    % --- Desired Sliding Variable Dynamics ---
    Sdot_des = -(c.c3 + c.c4*phi_as)*S - u_sta;

    % --- Equivalent Control Law Calculation ---
    other_plant_terms = p.k*theta_dot(k-1) + p.m*p.g*p.l*sin(theta(k-1));
    Sdot_deriv_terms = (c.c1 + c.c2*phi_a)*ed + c.c2*phi_a_dot*e;
    
    tau(k) = other_plant_terms - I * (Sdot_des - Sdot_deriv_terms);

    % --- Apply Actuator Limits ---
    tau(k) = max(-25, min(25, tau(k)));

    % --- Integrate Plant Dynamics ---
    theta_ddot   = (tau(k) - other_plant_terms) / I;
    theta_dot(k) = theta_dot(k-1) + dt*theta_ddot;
    theta(k)     = theta(k-1)     + dt*theta_dot(k); % Semi-implicit Euler
end
% Ensure first value of setpoint log is correct for plotting
setpoint_log_deg(1) = ref1_deg;

%% 4. PLOTTING RESULTS
figure('Name','PT-STSMC - Angular Position');
hold on; grid on;
plot(t, theta*180/pi, 'b', 'LineWidth', 2);
plot(t, setpoint_log_deg, 'k--', 'LineWidth', 2, 'DisplayName', 'Setpoint');
xlabel('Time [s]');
ylabel('\theta [deg]');
title('PT-STSMC: Angular Position Tracking');
legend('Response', 'Setpoint', 'Location', 'SouthEast');

figure('Name','PT-STSMC - Control Torque');
hold on; grid on;
plot(t, tau, 'r', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('\tau [Nm]');
title('PT-STSMC: Control Torque');

%% 5. PRESCRIBED-TIME HELPER FUNCTIONS
function out = phi_fun(t,t0,T,p)
    if t>=t0 && t<t0+T
        out = mu_dot_fun(t0,T,t,p)/(mu_fun(t0,T,t,p));
    else
        out = p/T;
    end
end

function out = mu_fun(t0,T,t,p)
    if t>=t0 && t<t0+T
        out = (T/(t0+T-t))^p;
    else
        out = 1;
    end
end

function out = phi_dot_fun(t,t0,T,p)
    if t>=t0 && t<t0+T
        out = (p/T^2)*(mu_fun(t0,T,t,p)^(2+(1/p)));
    else
        out = 0;
    end
end

function out = mu_dot_fun(t0,T,t,p)
    if t>=t0 && t<t0+T
        out = (p/T)*(mu_fun(t0,T,t,p)^(1+(1/p)));
    else
        out = 0;
    end
end