function [Kp, Ki] = it2flc_gains(E, DE, UP, UI)
% IT2FLC_GAINS  Interval type-2 fuzzy self-tuner of the FOPI gains
%
%   Inputs : E  = K_E  * df          (Eq. 37)
%            DE = K_DE * d(df)/dt    (Eq. 37)
%   Outputs: Kp, Ki (Eq. 38)
%
%   - 5 triangular IT2 MFs per input on [-1, 1]: NB NM Z PM PB (Fig. 11)
%       upper MF height 1, lower MF height 0.5, same support
%   - 5 output sets S M B VB VVB with centres 0, U/4, U/2, 3U/4, U (Fig. 12, Eq. 33)
%   - Product implication, centre-of-sets type reduction (Karnik-Mendel),
%     centroid defuzzification  y = (yl + yr)/2
%   - Rule base: Table 3a (Kp) and Table 3b (Ki)

E  = min(max(E , -1), 1);
DE = min(max(DE, -1), 1);

c   = [-1 -0.5 0 0.5 1];          % MF centres  NB NM Z PM PB
hw  = 0.5;                         % MF half-width
hL  = 0.5;                         % lower-MF height (footprint of uncertainty)

muU_E  = max(0, 1 - abs(E  - c)/hw);   muL_E  = hL*muU_E;
muU_DE = max(0, 1 - abs(DE - c)/hw);   muL_DE = hL*muU_DE;

% firing intervals (product t-norm), rows = df, cols = d(df)/dt
fU = muU_E.' * muU_DE;
fL = muL_E.' * muL_DE;

% output set index: 1=S 2=M 3=B 4=VB 5=VVB
RP = [5 5 4 4 3;      % df = NB
      5 4 4 3 2;      % df = NM
      4 3 3 3 2;      % df = Z
      4 3 2 2 1;      % df = PM
      3 2 2 1 1];     % df = PB
RI = [1 1 2 3 3;
      1 2 3 3 3;
      2 2 3 4 4;
      2 3 3 4 5;
      3 3 4 5 5];

cen = [0 0.25 0.5 0.75 1];        % normalised centres (Eq. 33)

Kp = UP * km_cos(fL(:), fU(:), cen(RP(:)));
Ki = UI * km_cos(fL(:), fU(:), cen(RI(:)));
end

% -------------------------------------------------------------------------
function y = km_cos(fl, fu, yc)
% Karnik-Mendel centre-of-sets type reduction (exhaustive switch-point search)
% (rules with zero firing add no weight, so no filtering is needed;
%  at least one rule always fires because the inputs are saturated to [-1,1])
[yc, o] = sort(yc(:));
fl = fl(o); fu = fu(o);

cu  = cumsum(fu.*yc);  cuw = cumsum(fu);
cl  = cumsum(fl.*yc);  clw = cumsum(fl);
tu  = cu(end); tuw = cuw(end);
tl  = cl(end); tlw = clw(end);

% left end point: upper firing for i<=k, lower firing for i>k
numL = [0; cu]  + (tl  - [0; cl]);
denL = [0; cuw] + (tlw - [0; clw]);
yl = min(numL ./ max(denL, eps));

% right end point: lower firing for i<=k, upper firing for i>k
numR = [0; cl]  + (tu  - [0; cu]);
denR = [0; clw] + (tuw - [0; cuw]);
yr = max(numR ./ max(denR, eps));

y = 0.5*(yl + yr);
end
