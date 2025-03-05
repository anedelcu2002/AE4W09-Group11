function [Blade] = scale_chord_and_twist_for_blade(Blade,lambda_old,lambda_new,R_new,B_old,B_new)
arguments
    Blade
    lambda_old
    lambda_new
    R_new
    B_old = 3
    B_new = 3
end

R_old = Blade.Radius(end);

% Calculate the scaling factors.
lambda_factor = lambda_new/lambda_old;
R_factor = R_new/R_old;
cl_factor = 1;  % TODO: We could add scaling based on the airfoil lift coefficient, but not needed for now.
B_factor = B_new / B_old;


% See 'System design and scaling' powerpoint, slide 10 and 28.
% Chord ~= R / (B * cl * lambda^2)
Blade.Chord = Blade.Chord * R_factor / (B_factor * cl_factor * lambda_factor^2);

% Twist ~= 1 / lambda
Blade.Twist = Blade.Twist / lambda_factor;

end
