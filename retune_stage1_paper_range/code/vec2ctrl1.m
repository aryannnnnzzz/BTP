function C = vec2ctrl1(x, Tf)
%VEC2CTRL1  Map the stage-1 design vector to a controller struct.
%   x = [Kp Ki Kd lam mu Kvi Dvi]
C = struct('Kp',x(1),'Ki',x(2),'Kd',x(3),'lam',x(4),'mu',x(5), ...
           'Kvi',x(6),'Dvi',x(7),'Tf',Tf);
end
