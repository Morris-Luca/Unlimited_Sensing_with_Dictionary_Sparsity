function [X, I] = OMP(Y, dico, mode, modeparam)
% OMP - Orthogonal Matching Pursuit Algorithm
% Performs sparse coding using Orthogonal Matching Pursuit (OMP) with 
% direct matrix-on-vector operations for any dictionary.
%
% INPUTS:
%   Y         - Input signal (d x 1 vector, where d is the signal dimension)
%   dico      - Dictionary (d x K matrix, where K is the number of atoms)
%   mode      - Mode of operation:
%               'a' - Adaptive: Add atoms as long as max_i |<d_i, r>| > th * ||r||_2
%               's' - Fixed sparsity: Use fixed sparsity level S (stored in modeparam)
%               'e' - Error-based: Add atoms until relative error is below threshold
%   modeparam - Parameter for the selected mode (e.g., sparsity level or error threshold)
%
% OUTPUTS:
%   X         - Sparse coefficient vector (K x 1)
%   I         - Support indices (indices of selected atoms)
%
% NOTES:
%   - Based on Lorenzo Granai's code.
%   - Adapted by Karin Schnass, 2024.

    % Check input arguments
    if nargin < 2
        disp('Missing input arguments');
        return;
    end

    % Get dimensions of dictionary and signal
    [d, K] = size(dico);
    [dd, N] = size(Y);
    if d ~= dd || N ~= 1
        disp('Mismatch between dictionary size and signal size, or multiple signals provided');
        return;
    end

    % Set default mode and parameters if not provided
    if nargin < 4
        mode = 'a';
        modeparam = 4; 
    end

    % Initialize parameters based on the selected mode
    if mode == 'a'
        th = sqrt(2 * (log(2 * K) + log(modeparam)) / d); % Threshold for adaptive mode
        S = floor(min(K, d) / 2); % Maximum sparsity level
        err = 0; % No error constraint
    elseif mode == 's'
        S = max(modeparam, 1); % Fixed sparsity level
        th = 0; % No threshold
        err = 0; % No error constraint
    elseif mode == 'e'
        S = floor(min(K, d) / 2); % Maximum sparsity level
        th = 0; % No threshold
        err = modeparam; % Error threshold
    end

    % Initialize output variables
    X = zeros(K, 1); % Sparse coefficient vector
    if mode == 's'
        I = zeros(S, 1); % Support indices for fixed sparsity
    else
        I = zeros(K, 1); % Support indices for other modes
    end

    % Precompute the Gram matrix of the dictionary
    fullGram = dico' * dico;

    % Extract the signal
    y = Y; % Input signal
    normy = norm(y); % Norm of the signal
    supp = []; % Support set (indices of selected atoms)
    coeff = []; % Coefficients for selected atoms
    GramInvNew = zeros(S); % Inverse Gram matrix (updated iteratively)

    % Initialize residual
    res = y;
    normres = normy;

    % Initialize iteration counter
    k = 0;

    % Compute initial inner products
    ip = dico' * res;
    [maxip, i] = max(abs(ip));

    % Main OMP loop
    while (maxip > th * normres) && (k < S) && (normres > err * normy)
        i = i(1); % Select the index of the maximum inner product
        newatom = dico(:, i); % Corresponding atom

        if k == 0
            % Initialize Gram matrix inverse for the first atom
            GramInvNew(1, 1) = 1 / fullGram(i, i);
            k = 1;
            GramInv = GramInvNew(1, 1);
        else
            % Update the inverse of the Gram matrix
            Q = GramInv * fullGram(supp, i);
            m = 1 / (fullGram(i, i) - fullGram(supp, i)' * Q);

            GramInvNew(1:k, 1:k) = GramInv + Q * Q' * m;
            GramInvNew(1:k, k+1) = -Q * m;
            GramInvNew(k+1, 1:k) = -Q' * m;
            GramInvNew(k+1, k+1) = m;

            k = k + 1;
            GramInv = GramInvNew(1:k, 1:k);
        end

        % Update support set and coefficients
        supp = [supp i];
        coeff = GramInv * (dico(:, supp)' * y);

        % Update residual
        res = y - dico(:, supp) * coeff;
        normres = norm(res);

        % Compute new inner products
        ip = dico' * res;
        [maxip, i] = max(abs(ip));
    end

    % Store coefficients and support indices
    X(supp) = coeff;
    if mode == 's'
        I(:) = supp;
    else
        I(supp) = 1;
    end
end