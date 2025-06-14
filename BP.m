function [coeff, support] = BP(signal_in, F, iterations)
% BP - Basis Pursuit algorithm adapted
%
% Solves the Basis Pursuit (BP) problem to find a solution `coeff` such that:
%   signal_in = F * coeff
% with the minimum l1 norm. The problem is solved using Linear Programming (LP)
% techniques via Interior-Point Methods.
%
% INPUT:
%   signal_in   - 1D input signal (n x 1 vector)
%   F           - Dictionary matrix (n x p matrix, where columns are dictionary atoms)
%   iterations  - Maximum number of iterations for the optimization solver
%
% OUTPUT:
%   coeff       - Coefficients vector (p x 1 vector)
%   support     - Indices of non-zero coefficients in `coeff`
%
% REQUIREMENTS:
%   - MATLAB Optimization Toolbox
%   - MATLAB Wavelet Toolbox (optional, depending on usage)
%
% REFERENCES:
%   Zhang, Y., "Solving Large-Scale Linear Programs by Interior-Point Methods 
%   Under the MATLAB Environment."
%
% NOTES:
%   - This implementation is adapted from BP.m (May 2003) and BPv4 (Nov 2004).
%   - Original author: Lorenzo Granai, EPFL, Nov 2004.
%   - Adapted by Karin Schnass, 2024.

    % Get the size of the input signal
    n = size(signal_in, 1);
    
    % Check if the input signal is a column vector
    if size(signal_in, 2) ~= 1
        error('ERROR in BP: Input signal must be a column vector (n x 1).');
    end

    % Check if the dictionary matrix F has the correct dimensions
    if size(F, 1) ~= n
        error('ERROR in BP: Dictionary matrix F must have the same number of rows as the input signal.');
    end

    % Get the dimensions of the dictionary matrix
    [d, p] = size(F);
    K = p - d; % Number of Dirac atoms in the dictionary

    % Construct the augmented matrix for the linear programming problem
    % A = [F, -F] to handle positive and negative coefficients
    A = [F, -F];

    % Set options for the linear programming solver
    % 'LargeScale' enables large-scale optimization
    % 'TolFun' sets the tolerance for the objective function
    % 'MaxIter' sets the maximum number of iterations
    % 'Display' controls the verbosity of the solver
    optMyBP = optimset('LargeScale', 'on', ...
                       'TolFun', 1.0e-2, ...
                       'MaxIter', iterations, ...
                       'Display', 'off');

    % Solve the linear programming problem using linprog
    % Minimize the l1 norm of the coefficients (sum of absolute values)
    % Subject to: A * x = signal_in
    [x, fval, exitflag] = linprog(ones(1, 2 * p), [], [], A, signal_in, zeros(2 * p, 1), [], optMyBP);

    % Handle the solver's exit flag
    if exitflag < 0
        % No convergence to a solution
        coeff = zeros(p, 1);
        support = [];
    elseif exitflag == 0
        % Maximum number of iterations exceeded
        coeff = zeros(p, 1);
        support = [];
    else
        % Successful optimization
        % Extract the coefficients from the solution vector
        coeff = x(1:p) - x(p+1:2*p); % Compute the difference between positive and negative parts
        coeff = coeff(1:K);          % No sparsity constraint applied

        % Find the support (indices of non-zero coefficients)
        support = find(coeff ~= 0);
    end
end