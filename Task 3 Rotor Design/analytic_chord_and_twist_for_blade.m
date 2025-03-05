function [Blade] = analytic_chord_and_twist_for_blade(Blade,Airfoil,lambda)

cl = zeros(size(Blade.Radius));
alpha_rad = zeros(size(Blade.Radius));

% Loop through the airfoils in this blade.
for i = 1:length(Blade.IFoil)
    i_airfoil = Blade.IFoil(i);
    alpha_array = Airfoil.Alpha{i_airfoil};
    Cl_array = Airfoil.Cl{i_airfoil};
    Cd_array = Airfoil.Cd{i_airfoil};
    l_over_d_array = Cl_array ./ Cd_array;
    [~, l_over_d_max_index] = max(l_over_d_array);

    cl_opt = Cl_array(l_over_d_max_index);
    alpha_opt = deg2rad(alpha_array(l_over_d_max_index));

    if all(Cl_array == 0)
        % This is a cylindrical airfoil, so definng an optimal angle of
        % attack makes no sense.
        alpha_opt = nan;
    end
       
    % Assign optimal cl and alpha to the blade locations that use this
    % airfoil.
    cl(Blade.NFoil == i_airfoil) = cl_opt;
    alpha_rad(Blade.NFoil == i_airfoil) = alpha_opt;
end

[twist_rad, chord] = calc_ideal_chord_and_twist(lambda, alpha_rad, cl, Blade.Radius);
Blade.Chord = chord;
Blade.Twist = rad2deg(twist_rad);
end


function [twist, chord] = calc_ideal_chord_and_twist(lambda, alpha, cl, r, B)
% Calculates the optimal twist and chord distribution at r (r can be a
% vector).
arguments
    lambda (1,1) double  % Design tip speed ratio
    alpha (1,:) double   % Design angle of attack at r.
    cl (1,:) double      % Airfoil lift coefficient at r.
    r (1,:) double       % Vector along blade length where alpha is defined.
    B (1,1) {mustBeInteger} = 3  % Number of blades.
end

% Get the radius of the blade.
R = r(end);

% Method from System design and scaling powerpoint, slide 10:
% Calculate ideal twist distribution.
twist = (2/3) ./ (lambda * r/R) - alpha;
twist = twist - twist(end);  % Convention that twist at the tip is 0.
% twist = rad2deg(twist);

% Calculate chord distribution.
chord_over_R = (16/9) * pi ./ (B * cl * lambda^2) .* (r/R).^(-1);
chord = chord_over_R * R;

end
