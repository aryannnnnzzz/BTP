function S = stage2_defaults()
%STAGE2_DEFAULTS  Initial parameters of the slow SOC-aware outer stage.
%   tune_stage2.m overwrites s1, s2, a_min, a_em, a_hi, Kr, kh_min, kh_max.

% SOC membership breakpoints (Low: SOC<=s1 fully, fades out by s2; High mirrored)
S.s1      = 0.15;
S.s2      = 0.35;
% |df| membership (Small below f_small, Large above f_big)
S.f_small = 0.005;
S.f_big   = 0.05;
% gain-scaling consequents
S.a_min   = 0.1;    % SOC at limit, frequency quiet: use very little of the ESS
S.a_em    = 0.5;    % SOC at limit, frequency emergency: partial support
S.a_hi    = 1.2;    % SOC on the "useful" side: allow extra support
% SOC recovery (only when frequency is quiet)
S.soc_ref  = 0.5;
S.soc_db   = 0.1;   % no recovery inside soc_ref +/- soc_db
S.Kr       = 0.2;   % p.u. power per unit SOC error
S.Prec_max = 0.03;  % recovery power cap (p.u.)
% ESS -> secondary-control handover gain (1/s)
S.kh_min  = 0.05;   % SOC healthy: slow handover (ESS helps longer)
S.kh_max  = 0.3;    % SOC near limits: fast handover to the thermal unit
% slow-stage timing
S.Ts2     = 0.1;    % scheduler sample time (s)
S.tau_a   = 2.0;    % smoothing of alpha and kh (s)
S.tau_r   = 30.0;   % smoothing of recovery power (s)
end
