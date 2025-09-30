%% Root Finding Methods for Financial Modeling
clear;
clc;

% Problem: Find the interest rate that gives a future value of Ugx15,000,000
% for an investment of Ugx10,000,000 over 5 years with monthly compounding
% Where: FV = Future Value, PV = Present Value, r = interest rate, n =
% years, t = period
% FV = PV * (1 + r/12)^(12*t)

PV = 10000000; FV_target = 15000000; t = 5; n = 12;

% Define the function
f = @(r) PV * (1 + r/n)^(n*t) - FV_target; %function
f_prime = @(r) PV * n*t * (1 + r/n)^(n*t - 1) / n; % derivative

%% Newton-Raphson Method
fprintf('=== NEWTON-RAPHSON METHOD ===\n');
tic;
r0 = 0.05; % Initial value (5%)
tolerance = 1e-8; %accuracy(dps)
max_iter = 100;
errors_newton = zeros(1,max_iter); %array for storing the error outputs

for i = 1:max_iter
    f_val = f(r0);
    f_deriv = f_prime(r0);
    r_new = r0 - f_val/f_deriv;
    
    error = abs(r_new - r0);
    errors_newton(i) = error;
    
    if error < tolerance
        break;
    end
    r0 = r_new;
    iter0 = i;
end

time_newton = toc;
fprintf('Root: r = %.6f (%.4f%%)\n', r0, r0*100);
fprintf('Iterations: %d, Time: %.6f seconds\n', iter0, time_newton);

%% Secant Method
fprintf('=== SECANT METHOD ===\n');
tic;
r0_secant = 0.04; % First initial guess
r1_secant = 0.06; % Second initial guess
tolerance = 1e-8; % accuracy (dps)
max_iter = 100;

r_secant = [r0_secant, r1_secant];
errors_secant = [];

for i = 1:max_iter
    f0 = f(r_secant(end-1));
    f1 = f(r_secant(end));
    
    if abs(f1 - f0) < eps
        break;
    end
    
    r_new = r_secant(i+1) - f1 * (r_secant(i+1) - r_secant(i)) / (f1 - f0);
    r_secant(i+2) = r_new;
    
    error = abs(r_secant(i+2) - r_secant(i+1));
    errors_secant(i:i) = error;
    
    if error < tolerance
        break;
    end
end

time_secant = toc;
fprintf('Root: r = %.6f (%.4f%%)\n', r_secant(end), r_secant(end)*100);
fprintf('Iterations: %d, Time: %.6f seconds\n', length(r_secant)-2, time_secant);

%% Visualization - Root Finding Comparison
figure;
% Function plot
r_range = linspace(0.01, 0.1, 100);
f_values = arrayfun(f, r_range);
plot(r_range, f_values, 'b-', 'LineWidth', 2);
hold on;
plot(r0, f(r0), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
plot(r_secant(end), f(r_secant(end)), 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'green');
xlabel('Interest Rate (r)');
ylabel('f(r)');
title('Function: f(r) = FV - FV_{target}');
legend('f(r)', 'Newton-Raphson Root', 'Secant Root', 'Location', 'best');
grid on;
saveas(gcf,'function_plot.png');

% Convergence comparison
figure;
semilogy(1:length(errors_newton), errors_newton, 'r-o', 'LineWidth', 2);
hold on;
semilogy(1:length(errors_secant), errors_secant, 'g-s', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Error (log scale)');
title('Convergence Comparison');
legend('Newton-Raphson', 'Secant', 'Location', 'best');
grid on;
saveas(gcf,'convergence_comparison.png');

% Verification
figure;
r_analytical = fzero(f, 0.05); % MATLAB's built-in function
methods = {'Newton-Raphson', 'Secant', 'fzero'};
roots = [r0, r_secant(end), r_analytical];
times = [time_newton, time_secant, 0.001]; % Approximate fzero time
bar(roots*100);
set(gca, 'XTickLabel', methods);
ylabel('Interest Rate (%)');
title('Root Comparison');
grid on;
saveas(gcf,'Verification.png');

% Add value labels on bars
for i = 1:length(roots)
    text(i, roots(i)*100 + 0.1, sprintf('%.4f%%', roots(i)*100), ...
        'HorizontalAlignment', 'center');
end

% Time comparison
figure;
bar(times*1000); % Convert to milliseconds
set(gca, 'XTickLabel', methods);
ylabel('Computation Time (ms)');
title('Computation Time Comparison');
grid on;
saveas(gcf,'Time_comparison.png');

% Add value labels on bars
for i = 1:length(times)
    text(i, times(i)*1000 + 0.1, sprintf('%.3f ms', times(i)*1000), ...
        'HorizontalAlignment', 'center');
end

% Final verification
figure;
FV_newton = PV * (1 + r0/n)^(n*t);
FV_secant = PV * (1 + r_secant(end)/n)^(n*t);
FV_analytical = PV * (1 + r_analytical/n)^(n*t);
bar([FV_newton, FV_secant, FV_analytical]);
hold on;
plot([0, 4], [FV_target, FV_target], 'r--', 'LineWidth', 2);
set(gca, 'XTickLabel', methods);
ylabel('Future Value ($)');
title('Future Value Verification');
legend('Calculated FV', 'Target FV', 'Location', 'best');
grid on;
saveas(gcf,'Final_verification.png');
sgtitle('Root Finding Methods: Compound Interest Rate Calculation', 'FontSize', 14, 'FontWeight', 'bold');


%% Real-world problem: Population growth with carrying capacity (Logistic Equation)
% dP/dt = r*P*(1 - P/K)
% Where: P = population, r = growth rate, K = carrying capacity

% Parameters
r = 0.1;        % Growth rate
K = 1000;       % Carrying capacity
P0 = 100;       % Initial population
t_span = [0, 50]; % Time span
h = 0.1;        % Step size

% Analytical solution
P_analytical = @(t) K ./ (1 + ((K - P0) / P0) * exp(-r * t));

%% Euler's Method
fprintf('=== EULER''S METHOD ===\n');
tic;
t_euler = t_span(1):h:t_span(2);
P_euler = zeros(size(t_euler));
P_euler(1) = P0;

for i = 1:length(t_euler)-1
    dPdt = r * P_euler(i) * (1 - P_euler(i)/K);
    P_euler(i+1) = P_euler(i) + h * dPdt;
end
time_euler = toc;

%% Runge-Kutta
fprintf('=== RUNGE-KUTTA  ===\n');
tic;
t_rk = t_span(1):h:t_span(2);
P_rk = zeros(size(t_rk));
P_rk(1) = P0;

for i = 1:length(t_rk)-1
    k1 = r * P_rk(i) * (1 - P_rk(i)/K);
    k2 = r * (P_rk(i) + 0.5*h*k1) * (1 - (P_rk(i) + 0.5*h*k1)/K);
    k3 = r * (P_rk(i) + 0.5*h*k2) * (1 - (P_rk(i) + 0.5*h*k2)/K);
    k4 = r * (P_rk(i) + h*k3) * (1 - (P_rk(i) + h*k3)/K);
    
    P_rk(i+1) = P_rk(i) + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
end
time_rk4 = toc;

% Error Analysis
P_analytical_values = P_analytical(t_euler);

error_euler = abs(P_euler - P_analytical_values);
error_rk = abs(P_rk(1:length(t_euler)) - P_analytical_values);

fprintf('=== ERROR ANALYSIS ===\n');
fprintf('Maximum Error - Euler: %.6f\n', max(error_euler));
fprintf('Maximum Error - RK4: %.6f\n', max(error_rk));
fprintf('Computation Time - Euler: %.6f seconds\n', time_euler);
fprintf('Computation Time - Runge kutta: %.6f seconds\n', time_rk4);

%% Plotting Results for Differential Equations

% Population growth comparison
figure;
t_dense = linspace(t_span(1), t_span(2), 1000);
plot(t_dense, P_analytical(t_dense), 'k-', 'LineWidth', 3, 'DisplayName', 'Analytical');
hold on;
plot(t_euler, P_euler, 'r--', 'LineWidth', 2, 'DisplayName', 'Euler');
plot(t_rk, P_rk, 'b:', 'LineWidth', 2, 'DisplayName', 'Runge kutta');
xlabel('Time');
ylabel('Population');
title('Population Growth: Logistic Equation');
legend('Location', 'southeast');
grid on;
saveas(gcf,'population_growth_comparison.png');

% Error comparison
figure;
semilogy(t_euler, error_euler, 'r--', 'LineWidth', 2, 'DisplayName', 'Euler Error');
hold on;
semilogy(t_euler, error_rk, 'b:', 'LineWidth', 2, 'DisplayName', 'Runge kutta Error');
xlabel('Time');
ylabel('Absolute Error (log scale)');
title('Error Comparison');
legend('Location', 'northeast');
grid on;
saveas(gcf,'error_comparison.png');

% Computation time
figure;
methods_de = {'Euler', 'Runge kutta'};
times_de = [time_euler, time_rk4] * 1000;
bar(times_de);
ylabel('Computation Time (ms)');
title('Computation Time');
set(gca, 'XTickLabel', methods_de);
grid on;
saveas(gcf,'computation_time.png');

% Phase portrait (dP/dt vs P)
figure;
P_range = linspace(0, 1200, 100);
dPdt_range = r * P_range .* (1 - P_range/K);
plot(P_range, dPdt_range, 'LineWidth', 2);
xlabel('Population (P)');
ylabel('dP/dt');
title('Phase Portrait');
grid on;
saveas(gcf,'phase_potrait.png');

% Final values comparison
figure;
final_values = [P_analytical(t_span(2)), P_euler(end), P_rk(end)];
bar(final_values);
ylabel('Final Population');
title('Final Population Values');
set(gca, 'XTickLabel', {'Analytical', 'Euler', 'Runge kutta'});
grid on;
saveas(gcf,'Final_comparison.png');

% Relative error over time
figure;
relative_error_euler = error_euler ./ P_analytical_values;
relative_error_rk = error_rk ./ P_analytical_values;
plot(t_euler, relative_error_euler * 100, 'r--', 'LineWidth', 2);
hold on;
plot(t_euler, relative_error_rk * 100, 'b:', 'LineWidth', 2);
xlabel('Time');
ylabel('Relative Error (%)');
title('Relative Error Over Time');
legend('Euler', 'Runge kutta', 'Location', 'northeast');
grid on;
saveas(gcf,'relative_error.png');
sgtitle('Differential Equation Solvers: Population Growth Model');