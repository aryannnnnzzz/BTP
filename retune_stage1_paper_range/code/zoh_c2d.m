function [Ad, Bd] = zoh_c2d(A, B, dt)
%ZOH_C2D  Exact zero-order-hold discretisation (no Control Toolbox needed).
n = size(A, 1);
if n == 0
    Ad = zeros(0); Bd = zeros(0, 1);
    return
end
Mx = expm([A, B; zeros(1, n + 1)] * dt);
Ad = Mx(1:n, 1:n);
Bd = Mx(1:n, n + 1);
end
