function P = mg_params()
%MG_PARAMS  Islanded-microgrid parameters (Sajadinia 2024, Table 2) plus the
%   ESS energy/SOC data that the base paper does not model.
%
%   GRC note: the paper lists "GRC = 20 %". Interpreted as 0.2 p.u./s this
%   script model reproduces Table 8 (without-VI and traditional-VI rows) to
%   within ~2 %. If your Step-1 Simulink model used another value, set it here.

% ---- microgrid (Table 2) ----
P.H     = 0.083;    % system inertia (p.u. MW s)
P.D     = 0.015;    % load damping (p.u. MW/Hz)
P.R     = 2.4;      % droop (Hz/p.u. MW)
P.Tg    = 0.1;      % governor time constant (s)
P.Tt    = 0.4;      % turbine time constant (s)
P.GRC   = 0.2;      % generation rate constraint (p.u./s)
P.VU    = 0.3;      % valve upper limit (p.u.)
P.VL    = -0.3;     % valve lower limit (p.u.)
P.beta  = 1.0;      % frequency bias (p.u. MW/Hz)
P.Ks    = 0.05;     % secondary integral gain
P.Twt   = 1.5;      % wind time constant (s)
P.Tpv   = 1.8;      % PV time constant (s)
P.Tinv  = 10;       % ESS / inverter time constant (s)
P.Pmax  = 0.3;      % ESS discharge limit (p.u.)
P.Pmin  = -0.3;     % ESS charge limit (p.u.)

% ---- ESS energy model (new in this project) ----
P.Sbase    = 20;    % MW (system base)
P.E_MWh    = 1.0;   % usable ESS energy (MWh); 5 MW / 1 MWh = 12-min battery
P.eta      = 0.95;  % one-way efficiency
P.soc_hmin = 0.05;  % BMS hard lower limit
P.soc_hmax = 0.95;  % BMS hard upper limit
P.soc_ramp = 0.01;  % BMS power taper width near the hard limits
P.soc_lo_ok = 0.10; % "healthy band" used in the SOC-violation metric
P.soc_hi_ok = 0.90;
end
