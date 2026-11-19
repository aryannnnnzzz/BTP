function S = vec2stage2(x, S)
%VEC2STAGE2  Map the stage-2 design vector onto a scheduler struct.
%   x = [s1 gap a_min a_em a_hi Kr kh_min dkh]   (s2 = s1+gap, kh_max = kh_min+dkh)
S.s1     = x(1);
S.s2     = min(x(1) + x(2), 0.5);
S.a_min  = x(3);
S.a_em   = max(x(4), x(3));
S.a_hi   = x(5);
S.Kr     = x(6);
S.kh_min = x(7);
S.kh_max = x(7) + x(8);
end
