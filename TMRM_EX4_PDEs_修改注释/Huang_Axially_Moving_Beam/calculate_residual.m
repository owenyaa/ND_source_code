function residual = calculate_residual(harmonic_coefficients)
    % Declare global variables for system properties
    global num_degrees_freedom num_harmonics time_data 
    global mass_matrix damping_matrix stiffness_matrix k12 k13 k22 k23 frequency force_amplitude1 force_amplitude2

    % Initialize displacement, velocity, and acceleration arrays
    displacement = zeros(num_degrees_freedom, length(time_data)); 
    velocity = zeros(num_degrees_freedom, length(time_data)); 
    acceleration = zeros(num_degrees_freedom, length(time_data));

    % Calculate displacements, velocities, and accelerations
    for dof_index = 1:num_degrees_freedom
        for harm_index = 1:num_harmonics
            harmonic_frequency = (2 * harm_index - 1) * frequency * time_data;
            displacement(dof_index, :) = displacement(dof_index, :) + harmonic_coefficients(harm_index, 2 * dof_index - 1) * cos(harmonic_frequency) + harmonic_coefficients(harm_index, 2 * dof_index) * sin(harmonic_frequency);
            velocity(dof_index, :) = velocity(dof_index, :) - frequency * (2 * harm_index - 1) * harmonic_coefficients(harm_index, 2 * dof_index - 1) * sin(harmonic_frequency) + frequency * (2 * harm_index - 1) * harmonic_coefficients(harm_index, 2 * dof_index) * cos(harmonic_frequency);
            acceleration(dof_index, :) = acceleration(dof_index, :) - (frequency * (2 * harm_index - 1))^2 * harmonic_coefficients(harm_index, 2 * dof_index - 1) * cos(harmonic_frequency) - (frequency * (2 * harm_index - 1))^2 * harmonic_coefficients(harm_index, 2 * dof_index) * sin(harmonic_frequency);
        end
    end

    % Compute system residual
    residual(1:num_degrees_freedom,:) = mass_matrix * acceleration + damping_matrix * velocity + stiffness_matrix * displacement + ...
                      [k12 * displacement(1, :) .* displacement(2, :).^2 + k13 * displacement(1, :).^3 - force_amplitude1 * cos(frequency * time_data); 
                       k22 * displacement(2, :) .* displacement(1, :).^2 + k23 * displacement(2, :).^3 - force_amplitude2 * cos(frequency * time_data)];

    % Calculate sensitivity of harmonic coefficients
    for sensitivity_index = 1:2 * num_harmonics * num_degrees_freedom
        sensitivity_params = zeros(2 * num_harmonics * num_degrees_freedom, 1);
        displacement_harm_sensitivity = zeros(num_degrees_freedom, length(time_data)); 
        velocity_harm_sensitivity = zeros(num_degrees_freedom, length(time_data)); 
        acceleration_harm_sensitivity = zeros(num_degrees_freedom, length(time_data));

        % Set current sensitivity coefficient to 1
        sensitivity_params(sensitivity_index, 1) = 1;
        reshaped_sensitivity = reshape(sensitivity_params, 2, num_harmonics * num_degrees_freedom);
        reshaped_sensitivity = reshaped_sensitivity';

        % Assign sensitivity parameters
        harmonic_sensitivity = reshaped_sensitivity(1:num_harmonics, 1:2);
        for dof_index = 1:num_degrees_freedom - 1
            harmonic_sensitivity = [harmonic_sensitivity, reshaped_sensitivity(dof_index * num_harmonics + 1:(dof_index + 1) * num_harmonics, 1:2)];
        end

        % Compute sensitivity displacements, velocities, and accelerations
        for dof_index = 1:num_degrees_freedom
            for harm_index = 1:num_harmonics
                harmonic_frequency = (2 * harm_index - 1) * frequency * time_data;
                displacement_harm_sensitivity(dof_index, :) = displacement_harm_sensitivity(dof_index, :) + harmonic_sensitivity(harm_index, 2 * dof_index - 1) * cos(harmonic_frequency) + harmonic_sensitivity(harm_index, 2 * dof_index) * sin(harmonic_frequency);
                velocity_harm_sensitivity(dof_index, :) = velocity_harm_sensitivity(dof_index, :) - frequency * (2 * harm_index - 1) * harmonic_sensitivity(harm_index, 2 * dof_index - 1) * sin(harmonic_frequency) + frequency * (2 * harm_index - 1) * harmonic_sensitivity(harm_index, 2 * dof_index) * cos(harmonic_frequency);
                acceleration_harm_sensitivity(dof_index, :) = acceleration_harm_sensitivity(dof_index, :) - (frequency * (2 * harm_index - 1))^2 * harmonic_sensitivity(harm_index, 2 * dof_index - 1) * cos(harmonic_frequency) - (frequency * (2 * harm_index - 1))^2 * harmonic_sensitivity(harm_index, 2 * dof_index) * sin(harmonic_frequency);
            end
        end

        % Compute residual for sensitivity coefficients
        residual(num_degrees_freedom * sensitivity_index + 1:num_degrees_freedom * (sensitivity_index + 1), :) = mass_matrix * acceleration_harm_sensitivity + damping_matrix * velocity_harm_sensitivity + stiffness_matrix * displacement_harm_sensitivity + ...
            [k12 * displacement(2, :).^2 .* displacement_harm_sensitivity(1, :) + 2 * k12 * displacement(2, :) .* displacement(1, :) .* displacement_harm_sensitivity(2, :) + 3 * k13 * displacement(1, :).^2 .* displacement_harm_sensitivity(1, :);
             k22 * displacement(1, :).^2 .* displacement_harm_sensitivity(2, :) + 2 * k22 * displacement(1, :) .* displacement(2, :) .* displacement_harm_sensitivity(1, :) + 3 * k23 * displacement(2, :).^2 .* displacement_harm_sensitivity(2, :)];
    end

    % Transpose residual for output
    residual = residual';
end







