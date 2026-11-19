function [best, bestFit, conv] = eo_optimize(fobj, lb, ub, o)
%EO_OPTIMIZE  Equilibrium Optimizer (Faramarzi et al., 2020), as in the base paper.
%   [best, bestFit, conv] = eo_optimize(fobj, lb, ub, o)
%   o.N (30)  o.MaxIt (50)  o.a1 (2)  o.a2 (1)  o.GP (0.5)  o.V (1)
%   o.useParallel (false)  o.verbose (true)  o.seed ([])
d = struct('N',30,'MaxIt',50,'a1',2,'a2',1,'GP',0.5,'V',1, ...
           'useParallel',false,'verbose',true,'seed',[]);
if nargin < 4, o = struct(); end
fn = fieldnames(d);
for i = 1:numel(fn), if ~isfield(o, fn{i}), o.(fn{i}) = d.(fn{i}); end, end
if ~isempty(o.seed), rng(o.seed); end

lb = lb(:)'; ub = ub(:)'; dim = numel(lb); N = o.N;
C      = lb + rand(N, dim).*(ub - lb);
Ceq    = zeros(4, dim);
Ceqfit = inf(4, 1);
C_old  = C;  fit_old = inf(N, 1);
conv   = zeros(o.MaxIt, 1);

for it = 1:o.MaxIt
    C   = min(max(C, lb), ub);                 % put exceeded variables on the bound
    fit = zeros(N, 1);
    if o.useParallel
        parfor i = 1:N, fit(i) = fobj(C(i,:)); end
    else
        for i = 1:N, fit(i) = fobj(C(i,:)); end
    end
    % memory saving: a particle never gets worse than its previous position
    worse = fit > fit_old;
    C(worse,:) = C_old(worse,:);
    fit(worse) = fit_old(worse);
    C_old = C;  fit_old = fit;

    % equilibrium pool: 4 best distinct candidates + their average
    allC = [Ceq; C];  allF = [Ceqfit; fit];
    [allF, ord] = sort(allF);  allC = allC(ord,:);
    keep = [true; diff(allF) > 1e-12];
    allF = allF(keep);  allC = allC(keep,:);
    mm = min(4, numel(allF));
    Ceq(1:mm,:) = allC(1:mm,:);  Ceqfit(1:mm) = allF(1:mm);
    pool = [Ceq; mean(Ceq, 1)];

    tt = (1 - it/o.MaxIt)^(o.a2*it/o.MaxIt);
    for i = 1:N
        lam = rand(1, dim);  r = rand(1, dim);
        ceq = pool(randi(5), :);
        F   = o.a1*sign(r - 0.5).*(exp(-lam.*tt) - 1);
        GCP = 0.5*rand*(rand >= o.GP);
        G   = GCP*(ceq - lam.*C(i,:)).*F;
        C(i,:) = ceq + (C(i,:) - ceq).*F + (G./(lam*o.V)).*(1 - F);
    end
    conv(it) = Ceqfit(1);
    if o.verbose
        fprintf('EO iter %3d/%d  best J = %.6g\n', it, o.MaxIt, Ceqfit(1));
    end
end
best = Ceq(1,:);  bestFit = Ceqfit(1);
end
