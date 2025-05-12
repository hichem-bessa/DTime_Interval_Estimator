function GG = MM_0(M, mu_, k)
%MM_0 Computes the polytopic matrix combination MM(ρ) at time step k.
%
% GG = MM_0(M, mu_, k)
%
% This function implements the MM(ρ) mapping as defined in:
%   "Finite-time estimation algorithms for LPV discrete-time systems with
%    application to output feedback stabilization"
%
% Inputs:
%   M   - Cell array of matrices {M_1, M_2, ..., M_{N+1}} where M{1} is the base term
%   mu_  - Cell array of scalar weighting functions {mu_1(k), mu_2(k), ..., mu_N(k)}
%   k   - Current time index or scalar parameter for evaluation
%
% Output:
%   GG  - Resulting matrix computed as:
%         GG = M{1} + sum_{i=1}^{N} M{i+1} * mu_{i}(k)

% Initialize with the first matrix
GG = M{1};

% Add weighted combinations from remaining matrices
for i = 2:length(mu_)+1
    GG = GG + M{i} * mu_{i-1}(k);
end

end
