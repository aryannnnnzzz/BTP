function sc = mg_scenarios(name, varargin)
%MG_SCENARIOS  Disturbance scenarios. Each returns
%   sc.T     : duration (s)
%   sc.dist  : @(t) -> n x 4  [dPresi dPind dPwind dPsolar]  (p.u.)
%   sc.soc0  : initial SOC
%   sc.E_MWh : [] (use mg_params) or an override
%
%   'I'          Fig. 13 random steps (values read off the figure - REPLACE
%                with the exact profile from your Step-1 model)
%   'II','III','IV'   Tables 7, 9, 11
%   'S5_lowSOC'  same-sign deficit build-up starting at SOC = 0.07
%   'S6_highSOC' mirrored surplus build-up starting at SOC = 0.93
%   'S7_long'    30-min random load/RES profile, SOC0 = 0.5 (drift, cycling)
%   'step'       single 0.1 p.u. industrial load step at t = 5 s, 600 s
%   Optional name/value: 'soc0', value  |  'E_MWh', value  |  'seed', value

ip = inputParser;
ip.addParameter('soc0', []);
ip.addParameter('E_MWh', []);
ip.addParameter('seed', 1);
ip.parse(varargin{:});

soc0 = 0.5;
switch upper(name)
    case 'I'
        T  = 100;
        Tb = (0:10:90)';
        % same profile as the Step-1 model (scenario_profile.m, id = 1);
        % the paper only gives Fig. 13, so both are read off that figure
        res = [0.021 0.001 0.003 0.022 0.012 0.027 0.028 0.012 0.026 0.021]';
        ind = [0.090 0.045 0.005 0.040 0.050 0.095 0.060 0.050 0.065 0.050]';
        wnd = [0.012 0.080 0.003 0.110 0.012 0.070 0.050 0.005 0.003 0.043]';
        pv  = [0.010 0.075 0.068 0.058 0.015 0.075 0.078 0.030 0.035 0.082]';
        M  = [res ind wnd pv];
    case 'II'
        T  = 50;  Tb = [0; 20; 40];
        M  = [0.05 0   0.12 0.08;
              0.05 0.1 0.12 0.08;
              0    0.1 0.12 0.08];
    case 'III'
        T  = 50;  Tb = [0; 20; 40];
        M  = [0.05 0.1 0.12 0   ;
              0.05 0.1 0.12 0.08;
              0.05 0.1 0    0.08];
    case 'IV'
        T  = 100; Tb = [0; 25; 45; 65; 85];
        M  = [0.05 0   0    0.08;
              0.05 0.1 0    0.08;
              0.05 0.1 0.12 0.08;
              0    0.1 0.12 0.08;
              0    0.1 0.12 0   ];
    case 'S5_LOWSOC'
        T  = 300; Tb = [0; 10; 60; 150; 200]; soc0 = 0.07;
        M  = [0    0   0.12 0.08;
              0.05 0   0.12 0.08;
              0.05 0.1 0.12 0.08;
              0.05 0.1 0    0.08;
              0.05 0.1 0    0   ];
    case 'S6_HIGHSOC'
        T  = 300; Tb = [0; 10; 60; 150; 200]; soc0 = 0.93;
        M  = [0.05 0.1 0    0   ;
              0    0.1 0    0   ;
              0    0   0    0   ;
              0    0   0.12 0   ;
              0    0   0.12 0.08];
    case 'S7_LONG'
        T  = 1800;
        rs = RandStream('mt19937ar', 'Seed', ip.Results.seed);
        Tb = 0; 
        while Tb(end) < T
            Tb(end+1,1) = Tb(end) + 5 + 25*rand(rs); %#ok<AGROW>
        end
        Tb = Tb(1:end-1);
        k  = numel(Tb);
        M  = [0.05*rand(rs,k,1), 0.10*rand(rs,k,1), ...
              0.12*rand(rs,k,1), 0.08*rand(rs,k,1)];
    case 'STEP'
        T  = 600; Tb = [0; 5];
        M  = [0 0 0 0; 0 0.1 0 0];
    otherwise
        error('mg_scenarios: unknown scenario "%s"', name);
end

sc.name  = name;
sc.T     = T;
sc.dist  = @(t) step_profile(Tb, M, t);
sc.soc0  = soc0;
sc.E_MWh = [];
if ~isempty(ip.Results.soc0),  sc.soc0  = ip.Results.soc0;  end
if ~isempty(ip.Results.E_MWh), sc.E_MWh = ip.Results.E_MWh; end
end

function D = step_profile(Tb, M, t)
idx = discretize(t(:), [Tb(:); inf]);
idx(isnan(idx)) = 1;
D = M(idx, :);
end
