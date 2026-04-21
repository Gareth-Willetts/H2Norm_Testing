function tfs = generateStableTFs()
    % Parameters
    degrees = 5:1:45;
    seed = 2;
    rng(seed); % Fixed seed for reproducibility
    
    executionTimes = zeros(size(degrees));
    tfs = cell(size(degrees)); % To store the resulting TF objects
    
    fprintf('Starting Transfer Function Generation...\n');
    
    for i = 1:length(degrees)
        n = degrees(i);
        
        % Start timing
        % tic;
        
        % Numerator: Normal distribution (mean 0, std 1)
        num = randn(1, n);
        
        % Denominator: Product of 2nd order factors (s^2 + a1*s + a0)
        % We need enough factors to reach degree n. 
        % If n is odd, we include one 1st order factor (s + a).
        denPoly = 1;
        remainingDegree = n;
        
        while remainingDegree > 0
            if remainingDegree >= 2
                % chi2rv with 2 d.o.f divided by 2
                % In MATLAB, chi2rnd(v) produces the distribution
                a_coeffs = chi2rnd(2, [1, 2]) / 2;
                % factor: s^2 + a1*s + a0
                denPoly = conv(denPoly, [1, a_coeffs(1), a_coeffs(2)]);
                remainingDegree = remainingDegree - 2;
            else
                % Handle odd degree with a first order factor (s + a0)
                a0 = chi2rnd(2) / 2;
                denPoly = conv(denPoly, [1, a0]);
                remainingDegree = remainingDegree - 1;
            end
        end
        
        % Create Transfer Function
        tfs{i} = tf(num, denPoly);
        
        % Record time
        % executionTimes(i) = toc;
        
        % Monitor first run and progress
        % if i == 1
        %     fprintf('First run (Degree %d) took: %.6f seconds\n', n, executionTimes(1));
        % else
        %     fprintf('Degree %d took: %.6f seconds\n', n, executionTimes(i));
        % end
    end
    
    % Plotting Execution Time
    % figure;
    % plot(degrees, executionTimes, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    % grid on;
    % title('Execution Time vs. Transfer Function Degree');
    % xlabel('Degree of Transfer Function');
    % ylabel('Execution Time (seconds)');
    
    % Verify stability of the last generated TF as a sample
    % fprintf('\nVerification: Poles of the degree %d system are in the LHP: %d\n', ...
    %     degrees(end), all(real(pole(tfs{end})) < 0));
end