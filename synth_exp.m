% Parameters
num_points = 10; % Number of points for sparsity and folding axes
num_averages = 1000; % Number of signal averages
success_tolerance = 0.05;
max_iterations = 2000; % Maximum iterations for BP

% Algorithm and dictionary types
algorithm_list = {'OMP', 'BP'};
dictionary_list = {'gdct', 'gabor_type',}; 
gap_frequency_list = [18, 24]; % For Gabor-type dictionary

% Loop over dictionary types
for dict_idx = 1:length(dictionary_list)
    dictionary_type = dictionary_list{dict_idx};
 
    % Dictionary-specific parameters
    if strcmp(dictionary_type, 'gdct')
        d = 1024;
        K = d / 2;
        folding_start = 0.3;
        folding_step = 0.02;
        sparsity_start = 10;
        sparsity_step = 10;
        loop_iterations = 1;
    elseif strcmp(dictionary_type, 'gabor_type')
        d = 1080;
        gap_time = 24;
        folding_step = 0.02;
        folding_start = 0.36;
        sparsity_start = 4;
        sparsity_step = 4;
        loop_iterations = 2;
    end

    % Calculate sparsity and folding axes
    sparsity_axis = (0:num_points-1) * sparsity_step + sparsity_start;
    folding_axis = (0:num_points-1) * folding_step + folding_start;
    
    % Loop over frequency gap for Gabor-type dictionaries 
    for gap_f_idx = 1:loop_iterations

        % Create dictionary
        if strcmp(dictionary_type, 'gdct')
            low_frequency_cutoff = 64;
    
            % Generate DCT dictionary
            dct_dictionary = idct(eye(d), 'type', 2);
            K = min(d - low_frequency_cutoff, K);
            dictionary = dct_dictionary(:, low_frequency_cutoff+1:low_frequency_cutoff + K);
    
            % Normalize the dictionary
            dictionary = dictionary * diag(1 ./ sqrt(sum(dictionary .* dictionary)));
        elseif strcmp(dictionary_type, 'gabor_type')
            gap_frequency = gap_frequency_list(gap_f_idx);
            time_vector = linspace(-d/2+1, d/2, d).';
            gaussian_window = exp(-pi * (time_vector / sqrt(d)).^2);
            gaussian_window = [gaussian_window(d/2+1:d); gaussian_window(1:d/2)];
    
            dictionary = make_gabor_type_dict(gaussian_window, gap_time, gap_frequency);
            K = size(dictionary, 2);
        end
    
        % Loop over algorithms
        for algo_idx = 1:length(algorithm_list)
            algorithm = algorithm_list{algo_idx};
    
            % Generate file path for saving results
            if strcmp(dictionary_type, 'gdct')
                numbering = 1 + (algo_idx - 1) * 3; % relevant for plotting order
                save_file = sprintf('%d_%s_%s_gauss_d_%d_K_%d_average_%d.mat', ...
                    numbering, algorithm, dictionary_type, d, K, num_averages);
            elseif strcmp(dictionary_type, 'gabor_type')
                numbering = algo_idx + 2 * gap_f_idx; % relevant for plotting order
                save_file = sprintf('%d_%s_%s_gauss_d_%d_average_%d_gapf_%d_gapt_%d.mat', ...
                    numbering, algorithm, dictionary_type, d, num_averages, gap_frequency, gap_time);
            end
    
            % Initialize result matrices
            error_matrix = ones(num_points, num_points);
            sparsity_matrix = ones(num_points, num_points);
            success_matrix = ones(num_points, num_points);
    
            % Loop over sparsity levels
            for sparsity_idx = 7:num_points
                current_sparsity = sparsity_axis(sparsity_idx);
    
                % Generate signals
                signals = [];
                support = zeros(current_sparsity, num_averages);
                for avg_idx = 1:num_averages
                    rng(avg_idx);
    
                    % Generate Gaussian-distributed sparse coefficients
                    coefficients = zeros(K, 1);
                    nonzero_indices = randperm(K, current_sparsity);
                    coefficients(nonzero_indices) = randn(current_sparsity, 1);
    
                    % Generate signal using the dictionary
                    signal = dictionary * coefficients;
    
                    % Append signal and support
                    signals = [signals, signal];
                    support(:, avg_idx) = nonzero_indices';
                end
    
                % Loop over folding heights
                for folding_idx = 1:num_points
                    folding_height = folding_axis(folding_idx);
    
                    % Perform folding operation
                    folding_operation = signals / (2 * folding_height) + 0.5;
                    cut_signals = floor(folding_operation);
                    sensed_signals = (folding_operation - cut_signals - 0.5) * 2 * folding_height;
    
                    % Calculate total folds
                    total_folds = zeros(1, num_averages);
                    for sig_idx = 1:num_averages
                        cut_signal = cut_signals(:, sig_idx);
                        total_folds(sig_idx) = sum(cut_signal(1:end-1) ~= cut_signal(2:end));
                    end
                    average_folds = mean(total_folds);
    
                    % Augment signals and dictionary for finite differences
                    sensed_signals_aug = [sensed_signals; sensed_signals(1, :)];
                    fd_sensed_sigs = sensed_signals_aug(2:d+1, :) - sensed_signals;
    
                    dictionary_aug = [dictionary; dictionary(1, :)];
                    finite_diff_dictionary = dictionary_aug(2:d+1, :) - dictionary;
                    norms_fd_dictionary = sqrt(sum(finite_diff_dictionary .* finite_diff_dictionary));
                    finite_diff_dictionary = finite_diff_dictionary ./ norms_fd_dictionary;
    
                    % Recovery dictionary
                    recovery_dictionary = [finite_diff_dictionary, eye(d)];
    
                    % Recovery error
                    error = zeros(1, num_averages);
                    for sig_idx = 1:num_averages
                        if strcmp(algorithm, 'OMP')
                            [X, ~] = OMP(fd_sensed_sigs(:, sig_idx), recovery_dictionary, 'e', 0.01);
                        elseif strcmp(algorithm, 'BP')
                            [X, ~] = BP(fd_sensed_sigs(:, sig_idx), recovery_dictionary, max_iterations);
                        end
    
                        % Reconstruct signals
                        reconstructed_signals = dictionary * diag(1 ./ norms_fd_dictionary) * X(1:K, :);
    
                        % Calculate error
                        signal_error  = norm(reconstructed_signals - signals(:, sig_idx), 'fro') / norm(signals(:, sig_idx), 'fro');
                        signal_error  = min(1, signal_error);
                        error(sig_idx) = signal_error;
                    end
    
                    % Update result matrices
                    error_matrix(folding_idx, sparsity_idx) = mean(error)
                    sparsity_matrix(folding_idx, sparsity_idx) = average_folds + current_sparsity;
                    success_matrix(folding_idx, sparsity_idx) = sum(error < success_tolerance) / num_averages;
    
                    % Save results
                    save(save_file, 'error_matrix', 'sparsity_matrix', 'success_matrix');
                end
            end

            % Save results
            save(save_file, 'error_matrix', 'sparsity_matrix', 'success_matrix', ...
                'num_points', 'num_averages', 'success_tolerance', ...
                'max_iterations', 'algorithm', 'dictionary_type', 'd', ...
                'sparsity_axis', 'folding_axis');

        end
    end
end
