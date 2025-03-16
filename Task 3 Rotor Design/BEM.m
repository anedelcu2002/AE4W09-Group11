clc
clear
close all

%% Airfoil data

DU_airfoil.dat = load("DU_97W300.txt");
DU_airfoil.alpha = DU_airfoil.dat(:,1);
DU_airfoil.Cl = DU_airfoil.dat(:,2);
DU_airfoil.Cd = DU_airfoil.dat(:,3);

NACAMid.dat = load("NACA_644421.txt");
NACAMid.alpha = NACAMid.dat(:,1);
NACAMid.Cl = NACAMid.dat(:,2);
NACAMid.Cd = NACAMid.dat(:,3);

NACAOut.dat = load("NACA_643218.txt");
NACAOut.alpha = NACAOut.dat(:,1);
NACAOut.Cl = NACAOut.dat(:,2);
NACAOut.Cd = NACAOut.dat(:,3);

AirfoilSet = [DU_airfoil NACAMid NACAOut];

%% Rotor Parameters

B = 3; % Number of Blades
R = 63; % Insert correct rotor radius
Pitch = 0; % NREL 5 MW has 0 set pitch but this could change

Twist = @(r) 14*(1-(r/R)); % Analytical equation for our turbine
Chord = @(r) 3*(1-(r/R)) + 1; % Analytical equation for our turbine

%% Operational Specs

lambda = 8; %Design tip speed ratio
U0 = 7.5; %U_infinity, average wind speed
rho = 1.225;

sigma_r = @(r) (B*Chord(r))/(2*pi*r); % Chord solidity

%% BEM
mu = linspace(0.2,1,100); %First entry depends on where cylinder ends
rRoot = mu(1)*R;
dr = mu(2)-mu(1); % deltar
mu = mu + 0.5*dr; % Shift radius location to be exactly in middle of segments
mu(end) = []; % Remove last entry since that is not on the blade

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
        if mu(j) < 0.4 % Inboard
            Airfoil = AirfoilSet(1);
        elseif mu(j) < 0.7 && mu(j) >= 0.4
            Airfoil = AirfoilSet(2);
        else
            Airfoil = AirfoilSet(3);
        end
        UR = U0*(1-a(i,j)); % Obtain axial velocity at rotor
        UTang = omega*r*(1+aprime(i,j)); % Obtain rotor-induced velocity
        Uapp = sqrt(UR^2 + UTang^2);

        phi(i,j) = atand(UR/UTang); % Inflow angle
        alpha(i,j) = phi(i,j) - Twist(r) - Pitch; % AoA

        [Cl,Cd] = interpolate_polars(Airfoil,alpha(i,j)); % Interpolation to find polars
        Cx = Cl*cosd(phi(i,j)) + Cd*sind(phi(i,j));
        Cy = Cl*sind(phi(i,j)) - Cd*cosd(phi(i,j));
        Ct(i,j) = ((Uapp^2)*Cx*Chord(r)*B)/((U0^2)*2*pi*r);
        
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

    aprime(i+1,j) = ((Uapp^2)*Chord(r)*Cy*B*R)/(8*pi*(r^2)*(U0^2)*(1-a(i+1,j))*lambda);
%% Tolerance check
    aDiff = abs(a(i+1,j) - a(i,j));
    aPrimeDiff = abs(aprime(i+1,j) - aprime(i,j));
    cond = (aDiff <= 1e-5 && aPrimeDiff <= 1e-5);
    i = i + 1;
    end
    CQ(1,j) = 4*aprime(i,j)*(1-a(i,j))*lambda*mu(j);
    CP(1,j) = 4*a(i,j)*(1-a(i,j))^2;
    CT(1,j) = 4*a(i,j)*(1-a(i,j));
end