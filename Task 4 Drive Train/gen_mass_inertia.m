function [M_gen, J_gen] = gen_mass_inertia(r_s, l_s, rho)
    % function that calculates the generator mass and inertia
    M_gen = pi*r_s*r_s*l_s*rho;
    J_gen = 0.5*M_gen*r_s*r_s; %Homogeneous cillinder approximation??????
end