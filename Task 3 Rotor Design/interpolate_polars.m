function [Cl_interp, Cd_interp] = interpolate_polars(DU_airfoil, alfa)

    if iscell(DU_airfoil.alpha)
        alpha = cell2mat(DU_airfoil.alpha);
        Cl = cell2mat(DU_airfoil.Cl);
        Cd = cell2mat(DU_airfoil.Cd);
    else
        alpha = DU_airfoil.alpha;
        Cl = DU_airfoil.Cl;
        Cd = DU_airfoil.Cd;
    end

Cl_interp = interp1(alpha, Cl, alfa, 'spline');
Cd_interp = interp1(alpha, Cd, alfa, 'spline');

end
