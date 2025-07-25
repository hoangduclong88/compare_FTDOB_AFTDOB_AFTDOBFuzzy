clear all;   clc;   close all;

% Simulation settings
dt = 0.0001; T = 3; N = T/dt;
t = 0:dt:T-dt;

% Reference signals
x1d_vals   = sin(2*pi*t);
dx1d_vals  = 2*pi*cos(2*pi*t);
ddx1d_vals = -4*pi^2*sin(2*pi*t);

% System model
theta1 = 5; theta2 = 8;
f = @(x2) -theta2*x2;
g = theta1;
d_func = @(ti) 8*sin(10*ti) + 4*cos(25*ti);

% Controller gains
% === FTDOB ===
r1 = 50; r2 = 25; c2 = 10;
b1F = 10; b2F = 5;

% === AFTDOB ===
gamma1_0 = 100; gamma2_0 = 40;
mu1 = 50; mu2 = 30; c1 = 10;
rho1 = 10; rho2 = 50;
b1A = 20; b2A = 15;

% === Fuzzy controller gains ===
b1AF = 30; b2AF = 20;

% Initialize variables
xF = [1; 0]; xA = [1; 0]; xAF = [1; 0];
xhatF = [0; 0]; xhatA = [0.5; 0.5]; xhatAF = [0.5; 0.5];

U_F = zeros(1,N); U_A = zeros(1,N); U_AF = zeros(1,N);
E_F = zeros(1,N); E_A = zeros(1,N); E_AF = zeros(1,N);
dtrue = zeros(1,N); dhatF = zeros(1,N); dhatA = zeros(1,N); dhatAF = zeros(1,N);

e_prev = 0;  % For fuzzy logic

% === Fuzzy Inference System ===
fis = mamfis('Name','AFTDOBFuzzy');
fis = addInput(fis, [-2 2], 'Name', 'Error');
fis = addInput(fis, [-500 500], 'Name', 'dError');
fis = addOutput(fis, [50 250], 'Name', 'gamma1');
fis = addOutput(fis, [20 120], 'Name', 'gamma2');

% MFs
fis = addMF(fis,'Error','gaussmf',[0.3 -1],'Name','N');
fis = addMF(fis,'Error','gaussmf',[0.3  0],'Name','Z');
fis = addMF(fis,'Error','gaussmf',[0.3  1],'Name','P');
fis = addMF(fis,'dError','gaussmf',[100 -300],'Name','DN');
fis = addMF(fis,'dError','gaussmf',[100 0],'Name','DZ');
fis = addMF(fis,'dError','gaussmf',[100 300],'Name','DP');

fis = addMF(fis, 'gamma1', 'trimf', [50 100 150], 'Name', 'L');
fis = addMF(fis, 'gamma1', 'trimf', [100 150 200], 'Name', 'M');
fis = addMF(fis, 'gamma1', 'trimf', [150 200 250], 'Name', 'H');

fis = addMF(fis, 'gamma2', 'trimf', [20 45 70], 'Name', 'L');
fis = addMF(fis, 'gamma2', 'trimf', [45 70 95], 'Name', 'M');
fis = addMF(fis, 'gamma2', 'trimf', [70 95 120], 'Name', 'H');

ruleList = [1 1 3 3 1 1; 1 2 3 2 1 1; 1 3 2 2 1 1; 2 2 2 2 1 1;
            2 1 2 3 1 1; 2 3 2 1 1 1; 3 1 2 2 1 1; 3 2 3 2 1 1; 3 3 3 3 1 1];
fis = addRule(fis, ruleList);

% === Simulation loop ===
for i = 1:N
    ti = t(i);
    x1d = x1d_vals(i); dx1d = dx1d_vals(i); ddx1d = ddx1d_vals(i);
    d = d_func(ti); dtrue(i) = d;

    %----------------- FTDOB -----------------
    eF = xF(2) - xhatF(2);
    d_hatF = r1 * eF + r2 * tanh(c2 * eF);
    dhatF(i) = d_hatF;
    dxhatF = [xhatF(2); f(xF(2)) + g*U_F(max(i-1,1)) + d_hatF];
    xhatF = xhatF + dt*dxhatF;
    z1 = xF(1) - x1d;
    alpha = dx1d - b1F*z1;
    z2 = xF(2) - alpha;
    alpha_dot = ddx1d - b1F*(xF(2) - dx1d);
    u = ( -d_hatF + alpha_dot + theta2*xF(2) - b2F*z2 - z1 ) / g;
    dx = [xF(2); g*u + f(xF(2)) + d];
    xF = xF + dt*dx;
    E_F(i) = z1; U_F(i) = u;

    xFdraw(1:2, i) = xF;

    r1draw(i) = r1;
    r2draw(i) = r2;

    %----------------- AFTDOB -----------------
    eA = xA(2) - xhatA(2);
    gamma1 = gamma1_0 + mu1 * min(abs(eA), rho1);
    gamma2 = gamma2_0 + mu2 * min(eA^2, rho2);
    d_hatA = gamma1*eA + gamma2*tanh(c1*eA);
    dhatA(i) = d_hatA;
    dxhatA = [xhatA(2); f(xA(2)) + g*U_A(max(i-1,1)) + d_hatA];
    xhatA = xhatA + dt*dxhatA;
    z1 = xA(1) - x1d;
    alpha = dx1d - b1A*z1;
    z2 = xA(2) - alpha;
    alpha_dot = ddx1d - b1A*(xA(2) - dx1d);
    u = ( -d_hatA + alpha_dot + theta2*xA(2) - b2A*z2 - z1 ) / g;
    dx = [xA(2); g*u + f(xA(2)) + d];
    xA = xA + dt*dx;
    E_A(i) = z1; U_A(i) = u;

    xAdraw(1:2, i) = xA;

    gamma1draw(i) = gamma1;
    gamma2draw(i) = gamma2;

    %----------------- AFTDOB-Fuzzy -----------------
    e = xAF(2) - xhatAF(2);
    de = (e - e_prev)/dt; e_prev = e;
    out = evalfis(fis, [e de]); g1 = out(1); g2 = out(2);
    d_hatAF = g1 * e + g2 * tanh(e * 10);
    dhatAF(i) = d_hatAF;
    dxhatAF = [xhatAF(2); f(xAF(2)) + g*U_AF(max(i-1,1)) + d_hatAF];
    xhatAF = xhatAF + dt*dxhatAF;
    z1 = xAF(1) - x1d;
    alpha = dx1d - b1AF*z1;
    z2 = xAF(2) - alpha;
    alpha_dot = ddx1d - b1AF*(xAF(2) - dx1d);
    u = ( -d_hatAF + alpha_dot + theta2*xAF(2) - b2AF*z2 - z1 ) / g;
    dx = [xAF(2); g*u + f(xAF(2)) + d];
    xAF = xAF + dt*dx;
    E_AF(i) = z1; U_AF(i) = u;

    xAFdraw(1:2, i) = xAF;

    g1draw(i) = g1;
    g2draw(i) = g2;

end

font = 12;
line = 1.3;
solution = 1600;

% === Plot Results ===
drawx1 = figure;
plot(t, x1d_vals, 'r', t, xFdraw(1,:), 'g', t, xAdraw(1,:), 'm', t, xAFdraw(1,:), 'b', 'LineWidth', line);
legend('$x_{1d}$', 'FTDOB', 'AFTDOB', 'AFTDOB-Fuzzy', 'interpreter', 'latex', 'FontSize', font);
xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', font);
grid on;

exportgraphics(drawx1, "x1.jpg", "Resolution", solution)

errorz1 = figure;
plot(t, E_F, 'g', t, E_A, 'm', t, E_AF, 'b', 'LineWidth', line);
legend('FTDOB', 'AFTDOB', 'AFTDOB-Fuzzy', 'interpreter', 'latex', 'FontSize', font); 
xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', font); 
grid on;

exportgraphics(errorz1, "errorz1.jpg", "Resolution", solution)

control = figure;
plot(t, U_F, 'g', t, U_A, 'm', t, U_AF, 'b',  'LineWidth', line);
legend('FTDOB', 'AFTDOB', 'AFTDOB-Fuzzy', 'interpreter', 'latex', 'FontSize', font); 
xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', font); 
grid on;

exportgraphics(control, "control.jpg", "Resolution", solution)

disturbance = figure;
plot(t, dtrue - dhatF, 'g', t, dtrue - dhatA, 'm', t, dtrue - dhatAF, 'b',   'LineWidth', 1.3);
legend('FTDOB', 'AFTDOB', 'AFTDOB-Fuzzy', 'interpreter', 'latex', 'FontSize', font);
xlabel('Time [s]', 'interpreter', 'latex', 'FontSize', font); 
grid on;

exportgraphics(disturbance, "disturbance.jpg", "Resolution", solution)

%% Draw Fuzzy:

% === Plot Membership Functions ===
MFinput = figure;
subplot(2,1,1); plotmf(fis, 'input', 1);
xlabel('Error', 'interpreter', 'latex', 'FontSize', font); 
ylabel('Degree of Membership', 'interpreter', 'latex', 'FontSize', font);
grid on;

% Get the current plot handles
h = findobj(gca, 'Type', 'Line');
% Increase the LineWidth for all membership function plots
set(h, 'LineWidth', line);  % Set desired line width

subplot(2,1,2); plotmf(fis, 'input', 2);
xlabel('dError', 'interpreter', 'latex', 'FontSize', font); 
ylabel('Degree of Membership', 'interpreter', 'latex', 'FontSize', font);
grid on;

% Get the current plot handles
h = findobj(gca, 'Type', 'Line');
% Increase the LineWidth for all membership function plots
set(h, 'LineWidth', line);  % Set desired line width

exportgraphics(MFinput, "MFinput.jpg", "Resolution", solution)

% === Plot Membership Functions of Fuzzy Outputs ===

% Plot for gamma1 (output 1)
MFoutput= figure;
subplot(2,1,1); plotmf(fis, 'output', 1);
xlabel('gamma1', 'interpreter', 'latex', 'FontSize', font); 
ylabel('Degree of Membership', 'interpreter', 'latex', 'FontSize', font);
grid on;

% Get the current plot handles
h = findobj(gca, 'Type', 'Line');
% Increase the LineWidth for all membership function plots
set(h, 'LineWidth', line);  % Set desired line width

subplot(2,1,2); plotmf(fis, 'output', 2);
xlabel('gamma2', 'interpreter', 'latex', 'FontSize', font); 
ylabel('Degree of Membership', 'interpreter', 'latex', 'FontSize', font);
grid on;

% Get the current plot handles
h = findobj(gca, 'Type', 'Line');
% Increase the LineWidth for all membership function plots
set(h, 'LineWidth', line);  % Set desired line width

exportgraphics(MFoutput, "MFoutput.jpg", "Resolution", solution)

% === Plot Fuzzy Rule Surface (for gamma1 and gamma2) ===
% Plot fuzzy rule surface for gamma1 (output 1)
figure;
gensurf(fis, [1 2], 1);  % Inputs 1 & 2 → Output 1
title('Fuzzy Surface for $\gamma_1$', 'interpreter', 'latex');

% Plot fuzzy rule surface for gamma2 (output 2)
figure;
gensurf(fis, [1 2], 2);  % Inputs 1 & 2 → Output 2
title('Fuzzy Surface for $\gamma_2$', 'interpreter', 'latex');

