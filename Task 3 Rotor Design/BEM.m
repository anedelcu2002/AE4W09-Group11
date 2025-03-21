

%% Airfoil data
Cylinder.dat = load("Cylinder.txt");
Cylinder.alpha = Cylinder.dat(:,1);
Cylinder.Cl = Cylinder.dat(:,2);
Cylinder.Cd = Cylinder.dat(:,3);

DU_airfoil.dat = load("DU_97W300.txt");
DU_airfoil.alpha = DU_airfoil.dat(:,1);
DU_airfoil.Cl = DU_airfoil.dat(:,2);
DU_airfoil.Cd = DU_airfoil.dat(:,3);

Blend1.dat = load("Blend DU_NACA.txt");
Blend1.alpha = Blend1.dat(:,1);
Blend1.Cl = Blend1.dat(:,2);
Blend1.Cd = Blend1.dat(:,3);

NACAMid.dat = load("NACA_644421.txt");
NACAMid.alpha = NACAMid.dat(:,1);
NACAMid.Cl = NACAMid.dat(:,2);
NACAMid.Cd = NACAMid.dat(:,3);

Blend2.dat = load("Blend NACA_NACA.txt");
Blend2.alpha = Blend2.dat(:,1);
Blend2.Cl = Blend2.dat(:,2);
Blend2.Cd = Blend2.dat(:,3);

NACAOut.dat = load("NACA_643218.txt");
NACAOut.alpha = NACAOut.dat(:,1);
NACAOut.Cl = NACAOut.dat(:,2);
NACAOut.Cd = NACAOut.dat(:,3);

AirfoilSet = [Cylinder DU_airfoil Blend1 NACAMid Blend2 NACAOut];

figure
hold on
plot(DU_airfoil.alpha,DU_airfoil.Cl,"DisplayName","DU 97-W-300")
plot(Blend1.alpha,Blend1.Cl,"DisplayName","Blend DU NACA");
plot(NACAMid.alpha,NACAMid.Cl,"DisplayName","NACA 644421");
plot(Blend2.alpha,Blend2.Cl,"DisplayName","Blend NACA NACA");
plot(NACAOut.alpha,NACAOut.Cl,"DisplayName","NACA 643218");
xlim([-20 20])
grid on
legend("Location",'best')
xlabel("\alpha");
ylabel("C_l")

%% Rotor Parameters

B = 3; % Number of Blades
R = 143/2; % Insert correct rotor radius
Pitch = 0; % NREL 5 MW has 0 set pitch but this could change

data = load("BulgAirChordTwist.mat")
Twist = data.Blade.Twist;
Chord = data.Blade.Chord;

%% Operational Specs

lambda = 8; %Design tip speed ratio
U0 = 7.5; %U_infinity, average wind speed
rho = 1.225;


%% BEM
NREL5MW = load("..\FASTTool\NREL5MW.mat");
BulgAir = load("BulgAir.mat").BulgAir;

r_n = NREL5MW.Blade.Radius / NREL5MW.Blade.Radius(end);
r_n_bulgAir = BulgAir.Blade.Radius / BulgAir.Blade.Radius(end);
BulgAir.Blade.Radius = interp1(r_n_bulgAir, BulgAir.Blade.Radius, r_n);
BulgAir.Blade.NFoil = interp1(r_n_bulgAir, BulgAir.Blade.NFoil, r_n, 'nearest');
mu = BulgAir.Blade.Radius / BulgAir.Blade.Radius(end);

rFace = zeros(1,length(mu)+1);
rFace(1) = 0;
rFace(end) = R;
for i = 2:length(mu)
    rFace(i) = 0.5*mu(i)*R + 0.5*mu(i-1)*R;
end
for i = 1:length(mu)
    dr(i) = rFace(i+1)-rFace(i);
end

omega = (lambda*U0)/R; %[Hz]
a = 0.3 * ones(1,length(mu)); % Initialize axial induction factor along blade span
aprime = 0.01 * ones(1,length(mu)); % Initialize tangential induction factor along blade span

for j = 1:length(mu)
    i = 1;
    cond = false;
    while cond == 0
        if i >= 2
        a(i,j) = 0.5*a(i,j) + 0.5*a(i-1,j);
        aprime(i,j) = 0.5*aprime(i,j) + 0.5*aprime(i-1,j);
        else

        end
        r =  R*mu(j); % What radial section are we looking at
        if mu(j) < 0.2 % Cylinder
            Airfoil = AirfoilSet(1);
        elseif mu(j) < 0.35 && mu(j) >= 0.2  % DU 97-W-300
            Airfoil = AirfoilSet(2);
        elseif mu(j) < 0.45 && mu(j) >= 0.35 % Blend 1
            Airfoil = AirfoilSet(3);
        elseif mu(j) < 0.6 && mu(j) >= 0.45  % NACA 644
            Airfoil = AirfoilSet(4);
        elseif mu(j) < 0.8 && mu(j) >= 0.6   % Blend 2
            Airfoil = AirfoilSet(5);
        else                                 % NACA 643
            Airfoil = AirfoilSet(6);
        end
        UR = U0*(1-a(i,j)); % Obtain axial velocity at rotor
        UTang = omega*r*(1+aprime(i,j)); % Obtain rotor-induced velocity
        Uapp = sqrt(UR^2 + UTang^2);

        phi(i,j) = atand(UR/UTang); % Inflow angle
        alpha(i,j) = phi(i,j) - Twist(j) - Pitch; % AoA

        [Cl,Cd] = interpolate_polars(Airfoil,alpha(i,j)); % Interpolation to find polars
        Cx = Cl*cosd(phi(i,j)) + Cd*sind(phi(i,j));
        Cy = Cl*sind(phi(i,j)) - Cd*cosd(phi(i,j));
        Ct(i,j) = ((Uapp^2)*Cx*Chord(j)*B)/((U0^2)*2*pi*r);
        
        a(i+1,j) = 0.5*(1-sqrt(1-Ct(i,j)));

%% Glauert Correction
        Ct1 = 1.816;
        if a(i+1,j) > 0.5 % Glauert Correction for induction factors above 0.5
            Ct(i,j) = 1.816 - 4*(sqrt(1.816)-1)*(1-a(i+1,j));
        else
            % Ct = 4*a(i+1,j)*(1-a(i+1,j));
        end

        Ct2 = 2*sqrt(Ct1) - Ct1;

        if Ct(i,j) >= Ct2
           a(i+1,j) = 1 + (Ct(i,j) - Ct1)/(4*sqrt(Ct1)-4);           
        else
           % a(i+1,j) = 0.5 - (sqrt(1-Ct)/2);
        end

    aprime(i+1,j) = ((Uapp^2)*Chord(j)*Cy*B*R)/(8*pi*(r^2)*(U0^2)*(1-a(i+1,j))*lambda);
%% Tolerance check
    aDiff = abs(a(i+1,j) - a(i,j));
    aPrimeDiff = abs(aprime(i+1,j) - aprime(i,j));
    cond = (aDiff <= 1e-5 && aPrimeDiff <= 1e-5);
    i = i + 1;
    end

%% Results
    FAxialSpan(1,j) = Cx*0.5*rho*(Uapp^2)*Chord(j); 
    FTangSpan(1,j) = Cy*0.5*rho*(Uapp^2)*Chord(j);
    CQ(1,j) = 4*aprime(i,j)*(1-a(i,j))*lambda*mu(j);
    CP(1,j) = 4*a(i,j)*(1-a(i,j))^2;
    CT(1,j) = 4*a(i,j)*(1-a(i,j));
end

%% Values needed for structural integration
    FAxial = FAxialSpan(1,:).*dr; 
    FTang = FTangSpan(1,:).*dr;

figure
hold on
plot(mu,FAxialSpan,"DisplayName","Axial Force")
plot(mu,FTangSpan,"DisplayName","Tangential Force")
title("Force in each blade section along span")
xlabel("Spanwise point \mu")
ylabel("Force per unit span [N/m]")
legend
grid on