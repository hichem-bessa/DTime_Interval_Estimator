function [mu_, A, B] = PolytopicModel(sys, borne_rho, A_func, B_func)
%POLYTOPICMODEL Constructs polytopic matrices A_i and B_i for LPV models.
%
% [mu_, A, B] = PolytopicModel(sys, borne_rho, A_func, B_func)
%
% Inputs:
%   sys        - Structure containing the rho functions (sys.rho{1}, sys.rho{2}, ...)
%   borne_rho  - Bounds for each rho parameter, size N x 2
%                Each row: [rho_max, rho_min]
%   A_func     - Function handle to compute A(rho)
%   B_func     - (Optional) Function handle to compute B(rho)
%
% Outputs:
%   mu_         - Cell array of weighting functions mu_i(k)
%   A          - Cell array of matrices A_i corresponding to vertices of the polytope
%   B          - Cell array of matrices B_i (if B_func is provided), otherwise empty

% Number of rho parameters
N = size(borne_rho, 1);
numVertices = 2^N;

% Initialize output cells
A = cell(1, numVertices);
mu_ = cell(1, numVertices);
if nargin == 4
    B = cell(1, numVertices);
else
    B = [];
end

% Iterate over each vertex of the hypercube
for i = 1:numVertices
    % Binary index to determine vertex selection
    binIndex = reverse(dec2bin(i - 1, N));
    
    % Construct vertex value v and reference beta for each rho dimension
    v = zeros(1, N);
    beta = zeros(1, N);
    
    for j = 1:N
        if binIndex(j) == '0'
            v(j) = borne_rho(j, 2);  % lower bound
            beta(j) = borne_rho(j, 1);  % upper reference
        else
            v(j) = borne_rho(j, 1);  % upper bound
            beta(j) = borne_rho(j, 2);  % lower reference
        end
    end
    
    % Define mu_i(k) as a product of normalized differences
    mu_{i} = @(k) prod(arrayfun(@(j) ...
        abs(sys.rho{j}(k) - beta(j)) / (borne_rho(j,2) - borne_rho(j,1)), ...
        1:N));
    
    % Evaluate A and optionally B at the vertex v
    A{i} = A_func(v);
    if nargin == 4
        B{i} = B_func(v);
    end
end
end
