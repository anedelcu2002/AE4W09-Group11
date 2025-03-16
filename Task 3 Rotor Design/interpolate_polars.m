function [Cl_interp, Cd_interp] = interpolate_polars(DU_airfoil, alfa)

alpha = DU_airfoil.alpha;
Cl = DU_airfoil.Cl;
Cd = DU_airfoil.Cd;

Cl_interp = interp1(alpha, Cl, alfa, 'spline');
Cd_interp = interp1(alpha, Cd, alfa, 'spline');

end
