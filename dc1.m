function fval = dc(t,P,Cc)

o1 = P(1);
o2 = P(2);
c1 = P(3);
c2 = P(4);

kaPos = 1500; %micrometer^12 micromole^-4 s^-1
kaNeg = 28.8; %s^-1
dc1dt = kaNeg*o1 - kaPos*(Cc^4)*c1;

% do2dt
kbPos = 1500; %micrometer^9 micromol^-3 s^-1
kbNeg = 385.9; %s^-1
do2dt = kbPos*(Cc^3)*o1 - kbNeg*o2;

% dc2dt
kcPos = 1.75; %s^-1
kcNeg = 0.1; %s^-1
dc2dt = kcPos*o1 - kcNeg*c2;

do1dt = -dc1dt -do2dt -dc2dt;
fval = [do1dt;do2dt;dc1dt;dc2dt];
end