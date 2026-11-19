function J = cost_stage1(x, ctx)
%COST_STAGE1  Fitness of a stage-1 candidate x = [Kp Ki Kd lam mu Kvi Dvi].
%   1) robustness gate: with a ctx.Tm_robust frequency-measurement lag the
%      frozen small-signal loop must be stable AND every oscillatory mode
%      must have damping ratio >= ctx.zeta_min (the paper's gains fail this)
%   2) normalised ITAE(df) + w_rocof*ITAE(RoCoF) + w_u*ISU over ctx.scen
C = vec2ctrl1(x, ctx.Tf);

if ctx.Tm_robust > 0
    [z, mr] = loop_margin(ctx.P, C, 1.0, ctx.opt.kh_fixed, ctx.Tm_robust, ctx.opt);
    if ~(mr < ctx.max_re) || z < ctx.zeta_min
        % graded penalty so EO is pulled back towards the feasible region
        J = 100 + 10*max(mr, 0) + 100*max(ctx.zeta_min - z, 0);
        return
    end
end
% optional second gate: larger lag, looser damping requirement
if isfield(ctx, 'Tm_robust2') && ctx.Tm_robust2 > 0
    [z, mr] = loop_margin(ctx.P, C, 1.0, ctx.opt.kh_fixed, ctx.Tm_robust2, ctx.opt);
    if ~(mr < ctx.max_re) || z < ctx.zeta_min2
        J = 50 + 10*max(mr, 0) + 100*max(ctx.zeta_min2 - z, 0);
        return
    end
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
