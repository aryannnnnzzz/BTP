function C = ctrl_presets(name)
%CTRL_PRESETS  Stage-1 controller structs.
%   Structure (all presets):  e = -df_meas
%     u1  = Kp*e + Ki*s^(-lam)*e + Kd*s^(mu)*e          (FO-PID, 0<lam<=1, 0<=mu<1)
%     Pc  = Kvi * d(u1)/dt + Dvi * u1                    (virtual rotor)
%   d/dt is implemented as s/(Tf*s+1).  Use Tf = 1e-3 with dt = 1e-3 to mimic
%   Simulink's du/dt block (needed to match the paper's peak values).
C = struct('Kp',0,'Ki',0,'Kd',0,'lam',1,'mu',0,'Kvi',0,'Dvi',0,'Tf',0.01);
switch lower(name)
    case 'none'
        % no VI
    case 'trad'                       % traditional VI (Table 4)
        C.Kp = 1; C.Kvi = 0.5;
    case 'pso_pi'                     % Table 4. NOTE: in this structure it does not
        C.Kp = 3.4957; C.Ki = 22.0991; C.Kvi = 0.5;   % reproduce Table 8; use your
                                                        % Step-1 implementation for it.
    case 'ga_pi'                      % Table 4
        C.Kp = 113.663; C.Ki = 122.14; C.Kvi = 3.71; C.Dvi = 4.794;
    case 'eo_fopi'                    % Table 4
        C.Kp = 109.195; C.Ki = 131.124; C.Kvi = 4.862; C.Dvi = 4.931; C.lam = 0.7192;
    case 'ostit2f'                    % base paper's proposed controller (Table 5)
        % IT2 fuzzy self-tuner of Kp, Ki around the EO-FOPI gains (Step-1
        % model, 'operating point' form); K_VI, D_VI, lambda from EO-FOPI
        C.Kp = 109.195; C.Ki = 131.124; C.Kvi = 4.862; C.Dvi = 4.931; C.lam = 0.7192;
        C.type = 'it2f';
        C.KE = 0.5327; C.KDE = 0.4937; C.UP = 153.8223; C.UI = 178.9173;
        C.Kp0 = 109.195; C.Ki0 = 131.124;
    case 'stage1_init'                % a sane starting point for stage-1 design
        C.Kp = 50; C.Ki = 0; C.Kd = 5; C.mu = 0.5; C.Kvi = 2; C.Dvi = 2;
    otherwise
        error('ctrl_presets: unknown preset "%s"', name);
end
end
