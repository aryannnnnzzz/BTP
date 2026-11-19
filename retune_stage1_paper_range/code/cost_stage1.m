function J = cost_stage1(x, ctx)
%COST_STAGE1  Fitness of a stage-1 candidate x = [Kp Ki Kd lam mu Kvi Dvi].
%   1) robustness gate: with a ctx.Tm_robust frequency-measurement lag the
%      frozen small-signal loop must be stable AND every oscillatory mode
%      must have damping ratio >= ctx.zeta_min (the paper's gains fail this)
%   2) normalised ITAE(df) + w_rocof*ITAE(RoCoF) + w_u*ISU over ctx.scen
C = vec2ctrl1(x, ctx.Tf);

% run 3: the gates are checked for every gain scaling alpha that stage 2
% may apply (ctx.alpha_set, default alpha = 1 only), so gain scaling stays usable
if isfield(ctx, 'alpha_set'), aset = ctx.alpha_set; else, aset = 1; end

% numerical-resolution gate: with the exact VI form, df feeds straight through
% to the ESS power with gain G = alpha*Kvi*(Kp + Kd*wh^mu)/Tinv (p.u./Hz).
% The fixed-step simulator is only accurate while dt*G/(2H) stays small, so
% a design that only "works" because of the time step is rejected.
if isfield(ctx, 'stiff_max') && ctx.Tf == 0
    G  = max(aset) * C.Kvi * (C.Kp + C.Kd*ctx.opt.wh^C.mu) / ctx.P.Tinv;
    st = ctx.opt.dt * G / (2*ctx.P.H);
    if st > ctx.stiff_max
        J = 200 + 100*(st - ctx.stiff_max);
        return
    end
end
gates = [ctx.Tm_robust, ctx.zeta_min];
if isfield(ctx, 'Tm_robust2') && ctx.Tm_robust2 > 0
    gates(end+1,:) = [ctx.Tm_robust2, ctx.zeta_min2];
end
viol = 0;  mrw = -inf;
for g = 1:size(gates, 1)
    if gates(g,1) <= 0, continue, end
    for a = aset
        [z, mr] = loop_margin(ctx.P, C, a, ctx.opt.kh_fixed, gates(g,1), ctx.opt);
        mrw = max(mrw, mr);
        if ~(mr < ctx.max_re) || z < gates(g,2)
            viol = viol + max(gates(g,2) - z, 0) + 0.01;
        end
    end
end
if viol > 0
    % graded penalty so EO is pulled back towards the feasible region
    J = 100 + 10*max(mrw, 0) + 100*viol;
    return
end

J = 0;
for s = 1:numel(ctx.scen)
    out = sim_mg(ctx.P, C, [], ctx.scen{s}, ctx.opt);
    if out.unstable, J = 1e4; return, end
    r = ctx.ref(s);
    J = J + out.J.ITAE/r.ITAE ...
          + ctx.w_rocof*out.J.ITAE_rocof/r.ITAE_rocof ...
          + ctx.w_u*out.J.ISU/r.ISU;
end
J = J / numel(ctx.scen);
end
