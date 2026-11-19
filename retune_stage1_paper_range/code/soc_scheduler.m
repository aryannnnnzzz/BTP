function [a_dis, a_ch, kh, prec] = soc_scheduler(soc, fabs, S)
%SOC_SCHEDULER  Stage 2: zero-order Sugeno fuzzy scheduler (no toolbox needed).
%   Inputs : soc  - state of charge (0..1)
%            fabs - low-pass filtered |df| (Hz)
%   Outputs: a_dis - scaling of K_VI, D_VI when the ESS must DISCHARGE (Pc>0)
%            a_ch  - scaling of K_VI, D_VI when the ESS must CHARGE    (Pc<0)
%            kh    - ESS -> secondary-control handover gain (1/s)
%            prec  - SOC-recovery power request (p.u., +ve = discharge)
%   Codegen-friendly: can be pasted into a Simulink MATLAB Function block.

% --- SOC memberships: Low / Mid / High (s2 <= 0.5 assumed) ---
w   = max(S.s2 - S.s1, 1e-6);
muL = min(max((S.s2 - soc) / w, 0), 1);
muH = min(max((soc - (1 - S.s2)) / w, 0), 1);
muM = max(1 - muL - muH, 0);

% --- |df| memberships: Small / Large ---
nL = min(max((fabs - S.f_small) / (S.f_big - S.f_small), 0), 1);
nS = 1 - nL;

% --- rule base ---
% discharge: SOC Low -> weak (a_min if quiet, a_em if emergency); Mid -> 1; High -> a_hi
% charge   : mirrored (SOC High is the "weak" side)
weak  = nS * S.a_min + nL * S.a_em;
a_dis = muL * weak + muM * 1 + muH * S.a_hi;
a_ch  = muH * weak + muM * 1 + muL * S.a_hi;

% handover: slow when SOC is healthy, fast near the limits
kh = S.kh_min + (S.kh_max - S.kh_min) * (1 - muM);

% SOC recovery: only outside the dead-band and only when frequency is quiet
dev  = soc - S.soc_ref;
ex   = sign(dev) * max(abs(dev) - S.soc_db, 0);
prec = min(max(S.Kr * ex, -S.Prec_max), S.Prec_max) * nS;
end
