clear; close all;

% Initialize initial condition for the differential equation
initial_condition = 0;

% Define the time data for the simulation
time_data = 0:0.01:80;

% Solve the differential equation using 'sub_liao13' as the function
[time_values, solution] = ode45('sub_liao13', time_data, initial_condition);

% Plot the numerical solution of the differential equation
figure;
plot(time_values, solution(:, 1), 'r.', 'MarkerSize', 15);

% Add legend and format the plot
legend_handle = legend('$$x$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
