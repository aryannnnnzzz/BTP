function opt = sim_options()
%SIM_OPTIONS  Default numerical options for sim_mg.
opt.dt       = 1e-3;   % fixed step (s). Use 2e-3 for tuning, 1e-3 for final results
opt.Tm       = 0;      % frequency-measurement lag time constant (s), 0 = ideal
opt.Td       = 0;      % pure measurement/communication delay (s)
opt.wb       = 1e-3;   % Oustaloup lower band edge (rad/s)
opt.wh       = 1e3;    % Oustaloup upper band edge (rad/s)
opt.N        = 4;      % Oustaloup order (2N+1 sections)
opt.kh_fixed = 0;      % fixed ESS->secondary handover gain when stage 2 is off
opt.decim    = 10;     % store every decim-th sample (metrics use every sample)
end
