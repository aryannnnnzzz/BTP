function [zeta, maxre, ev] = loop_margin(P, C, alpha, kh, Tm, opt)
%LOOP_MARGIN  Stability margin of the frozen small-signal loop.
%   zeta  : smallest damping ratio of the oscillatory modes (|Im| > 0.5 rad/s)
%   maxre : largest real part of all eigenvalues
%   A loop can have maxre < 0 and still be almost undamped (zeta ~ 0), which
%   rings for a long time and turns into a limit cycle once the ESS saturates.
%   Tuning therefore requires BOTH maxre < 0 and zeta >= zeta_min.
if nargin < 6 || isempty(opt), opt = sim_options(); end
ev    = mg_linearize(P, C, alpha, kh, Tm, opt);
maxre = max(real(ev));
osc   = ev(abs(imag(ev)) > 0.5);
if isempty(osc)
    zeta = 1;
else
    zeta = min(-real(osc) ./ abs(osc));
end
end
