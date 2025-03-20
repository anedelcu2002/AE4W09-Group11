%TASK 2 DAT
P_r = 3.5*10^6; %W
Omega_LSS_max = 1.0327; %rad/s
Omega_LSS_min = 0.224; %rad/s

%CONVERTER DATA
f_min = 4; %Hz
f_max = 50; %Hz

%GENERATOR DATA
n = 4;
F_d = 20; %kN/m^2
r_s = 2; %m

[Omega_HSS_min_allowed, Omega_HSS_max_allowed] = Omega_HSS_allowed_range(n, f_min, f_max); %rad/s
f_HSS_min_allowed = Omega_HSS_min_allowed/(2*pi);
f_HSS_max_allowed = Omega_HSS_max_allowed/(2*pi);

r_gbx = gbx_ratio_total(Omega_LSS_max, Omega_LSS_min, Omega_HSS_max_allowed, Omega_HSS_min_allowed, 150)
stage_types = [true,false,false];
r_list = gbx_design(r_gbx, 3, stage_types);
r_total = prod(r_list)

[Omega_HSS_min, Omega_HSS_max] = Omega_HSS_used_range(r_total, Omega_LSS_min, Omega_LSS_max);
f_HSS_min = Omega_HSS_min/(2*pi);
f_HSS_max = Omega_HSS_max/(2*pi);

eff_gbx = gbx_eff(r_list,stage_types)

V_gen = gen_sizing(F_d, P_r, Omega_HSS_max);
l_s = gen_length(V_gen, r_s);

[M_gen, J_gen] = gen_mass_inertia(r_s, l_s, 1000);


