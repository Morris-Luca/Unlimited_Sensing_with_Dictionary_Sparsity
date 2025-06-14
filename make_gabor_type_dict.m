function gabor_dict = make_gabor_type_dict(atom, gapt, gapf)
% MAKE_GABOR_TYPE_DICT - Generates a Gabor frame (dictionary) from a given atom and lattice.
%
% Input:
%   atom         - Input vector (Gabor atom)
%   gapt, gapf   - Lattice constants (must be integer divisors of n)
%
% Output:
%   gabor_dict   - Gabor frame matrix (size (l*p) x n), where l*p is the
%                  number of points in the Time-Frequency plane (l=n/gapt, p=n/gapf)
%
% Usage:
%   make_gabor_type_dict(atom, gapt, gapf);
%   make_gabor_type_dict(atom, gap);  % In this case, gapt = gapf = gap
%
% Author(s): H.G. Feichtinger, 09/1992
% Copyright: (c) NUHAG, Dept. Math., University of Vienna, Austria
%            http://nuhag.mat.univie.ac.at/
%            Permission is granted to modify and re-distribute this
%            code in any manner as long as this notice is preserved.
%            All standard disclaimers apply.

% Display help if no input arguments are provided
if nargin == 0
    help make_gabor_type_dict;
    return;
end

% Ensure atom is a row vector
atom = atom(:).';

% If only one gap is provided, set both gapt and gapf to the same value
if (nargin == 2) || (nargin == 1)
    gapf = gapt;
end

% Initialize parameters
n = length(atom);          % Length of the atom
ht = n / gapt;             % Number of time shifts
yp = 1:gapf:n;             % Frequency indices
xp = 1:gapt:n;             % Time indices
gabor_dict = [];           % Initialize Gabor frame matrix
aa = [atom, atom];         % Periodic extension of the atom
bas = -n/2+1:n/2;          % Frequency base indices

% Generate modulation matrix
em = exp(2 * pi * 1i * (bas - 1) / n * gapf); % Complex exponential for modulation
EM = [ones(1, n)];                            % Initialize modulation matrix

% Determine the number of frequency shifts (hf)
if gapf == 18
    hf = 10;
elseif gapf == 24
    hf = 7;
end

% Define the range of frequency shifts
j_values = [1:hf, (n/gapf - hf):(n/gapf - 1)];

% Construct the modulation matrix for all frequency shifts
for j = j_values
    EM = [EM; em.^j];
end

% Generate the Gabor frame
for jj = 1:ht
    % Circularly shift the atom
    rota = aa((n + 2 - xp(jj)):((2 * n) + 1 - xp(jj)));
    
    % Modulate the shifted atom
    modulated = ((ones(size(EM, 1), 1) * rota) .* EM);
    
    % Combine real and imaginary parts
    gabor_dict = [gabor_dict; real(modulated) + imag(modulated)];
end

% Transpose the Gabor frame matrix
gabor_dict = gabor_dict.';

% Normalize the Gabor frame
norms_gabor_dict = sqrt(sum(gabor_dict .* gabor_dict));
gabor_dict = gabor_dict ./ norms_gabor_dict;