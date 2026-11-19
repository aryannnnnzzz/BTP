function [A, B, Cc, Dd] = oustaloup_ss(alpha, wb, wh, N)
%OUSTALOUP_SS  State-space Oustaloup approximation of s^alpha on [wb, wh].
%   Cascade of 2N+1 first-order sections (s+z_k)/(s+p_k) with gain wh^alpha.
%   Valid for -1 < alpha < 1. alpha == 0 returns an empty state (gain 1).
if nargin < 2, wb = 1e-3; end
if nargin < 3, wh = 1e3;  end
if nargin < 4, N  = 4;    end
if abs(alpha) < 1e-9
    A = zeros(0); B = zeros(0,1); Cc = zeros(1,0); Dd = 1;
    return
end
assert(abs(alpha) < 1, 'oustaloup_ss: |alpha| must be < 1');
k  = -N:N;
M  = 2*N + 1;
z  = wb * (wh/wb).^((k + N + 0.5*(1 - alpha)) / M);
p  = wb * (wh/wb).^((k + N + 0.5*(1 + alpha)) / M);
K  = wh^alpha;
A  = diag(-p);
for i = 2:M
    A(i, 1:i-1) = z(1:i-1) - p(1:i-1);
end
B  = ones(M, 1);
Cc = K * (z - p);
Dd = K;
end
