function J = cost_stage2(x, ctx)
%COST_STAGE2  Fitness of a stage-2 candidate x = [s1 gap a_min a_em a_hi Kr kh_min dkh]
%   with stage 1 frozen (time-scale separation).
S = vec2stage2(x, ctx.S0);

% robustness gate over the frozen-gain envelope
if ctx.Tm_robust > 0
    for a = [S.a_min, 1, S.a_hi]
        for kh = [S.kh_min, S.kh_max]
            [z, mr] = loop_margin(ctx.P, ctx.C1, a, kh, ctx.Tm_robust, ctx.opt);
            if ~(mr < ctx.max_re) || z < ctx.zeta_min, J = 100; return, end
        end
    end
end

J = 0;
for s = 1:numel(ctx.scen)
    out = sim_mg(ctx.P, ctx.C1, S, ctx.scen{s}, ctx.opt);
    if out.unstable, J = 1e4; return, end
    r = ctx.ref(s);
    Js = out.J.ITAE/r.ITAE ...
       + ctx.w_thr*out.J.throughput_pus/r.throughput_pus ...
       + ctx.w_soc*out.J.soc_viol;
    % hitting the BMS limit = the ESS drops out = hard failure
    if out.J.socmin <= ctx.P.soc_hmin + 0.005 || out.J.socmax >= ctx.P.soc_hmax - 0.005
        Js = Js + ctx.w_bms;
    end
    J = J + Js;
end
J = J / numel(ctx.scen);
end
