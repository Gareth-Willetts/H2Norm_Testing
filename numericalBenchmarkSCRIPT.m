% =========================================================================
% --- UPDATED BENCHMARK SCRIPT ---
% =========================================================================
% tfs: Cell array of transfer functions from generateStableTFs()

% Get degrees from the generated cell array
num_tfs = length(tfs);
degrees = zeros(1, num_tfs);

% Number of repetitions for each system to get a reliable average
NUM_REPETITIONS = 1000; 
custom_times = zeros(1, num_tfs);
matlab_times = zeros(1, num_tfs);

fprintf('Starting H2 Norm Algorithm Benchmark...\n');
fprintf('Testing %d repetitions per system.\n', NUM_REPETITIONS);
fprintf('--------------------------------------------------\n');

for i = 1:num_tfs
    tf_sys = tfs{i};
    
    % Extract coefficients for custom function
    % tf.num/den returns 1x1 cells for SISO systems
    cn = tf_sys.num{1}; 
    cd = tf_sys.den{1};
    if length(cn) - length(cd) == 0
        cn = cn(2:end);
    end
    
    % Record degree for plotting (degree of denominator)
    degrees(i) = length(cd) - 1;

    WARM_UP_REPS = 5;
    for k = 1:WARM_UP_REPS
        temp_custom = numerical_H2_norm(cd,cn);
        temp_matlab = norm(tf_sys,2)^2;
    end
    
    % --- Measure Custom Algorithm Speed ---
    tic;
    for j = 1:NUM_REPETITIONS
        H2ns_custom = numerical_H2_norm(cd, cn);
        % H2ns_custom
    end
    custom_times(i) = toc / NUM_REPETITIONS;
    
    % --- Measure MATLAB Built-in Speed ---
    tic;
    for j = 1:NUM_REPETITIONS
        % Note: norm(sys, 2) returns the H2 norm, so we square it
        H2ns_matlab = norm(tf_sys, 2)^2;
    end
    matlab_times(i) = toc / NUM_REPETITIONS;
    
    % --- Validation ---
    % Using a relative tolerance check
    if abs(H2ns_custom - H2ns_matlab) > (0.001 * H2ns_matlab)
        fprintf('WARNING: Accuracy mismatch at Degree %d!\n', degrees(i));
        % H2ns_custom
        % H2ns_matlab
    end
    
    fprintf('Degree n=%-3d | Custom Avg Time: %8.6f s | MATLAB Avg Time: %8.6f s\n', ...
            degrees(i), custom_times(i), matlab_times(i));
end

fprintf('--------------------------------------------------\n');
fprintf('Benchmark Complete.\n');

% --- Plotting Results ---
figure;
plot(degrees, custom_times * 1e3, '-o', 'LineWidth', 1.5, 'DisplayName', 'Numerical H2 Norm Algorithm');
hold on;
plot(degrees, matlab_times * 1e3, '-x', 'LineWidth', 1.5, 'DisplayName', 'MATLAB norm Function');
hold off;
title('H2 Norm Squared Calculation Time vs. Transfer Function Degree');
xlabel('Transfer Function Degree (n)');
ylabel('Average Execution Time (ms)');
legend('Location', 'northwest');
grid on;
set(gca, 'FontSize', 16)

mdl_custom = fitlm(degrees, custom_times * 1e3);
mdl_matlab = fitlm(degrees, matlab_times * 1e3);
mdl_matlab_quadratic = fitlm(degrees, matlab_times * 1e3, 'quadratic');
figure;
plot(mdl_custom);
figure;
plot(mdl_matlab);
figure;
plot(mdl_matlab_quadratic);
mdl_custom
mdl_matlab
mdl_matlab_quadratic

figure;
subplot(1,2,1); plot(mdl_custom); title("Numerical H2 Norm Algorithm");
xlabel("Transfer Function Degree (n)"); ylabel("Average Execution Time (ms)");
grid on; legend('Location', 'northwest'); set(gca, 'FontSize', 16)
subplot(1,2,2); plot(mdl_matlab_quadratic); title("MATLAB norm Function");
xlabel("Transfer Function Degree (n)"); ylabel("Average Execution Time (ms)")
grid on; legend('Location', 'northwest');
set(gca, 'FontSize', 16)


% pCustom = polyfit(degrees, custom_times * 1e3, 1);
% pMatlab = polyfit(degrees, matlab_times * 1e3, 2);
% figure;
% plot(degrees, polyval(pCustom, degrees));
% hold on; 
% plot(degrees, custom_times * 1e3, '-o', 'LineWidth', 1.5, 'DisplayName', 'Numerical H2 Norm Algorithm');
% figure;
% plot(degrees, polyval(pMatlab, degrees));
% hold on;
% plot(degrees, matlab_times * 1e3, '-x', 'LineWidth', 1.5, 'DisplayName', 'MATLAB norm(tf, 2)^2');

%pCustom
%pMatlab

% Display Summary
ratio = custom_times ./ matlab_times;
[max_ratio, idx] = max(ratio);
fprintf('\nSummary: Custom algorithm is %.2fx slower than MATLAB at degree %d.\n', ...
        max_ratio, degrees(idx));

% figure;
% plot(degrees, custom_times*1e3)