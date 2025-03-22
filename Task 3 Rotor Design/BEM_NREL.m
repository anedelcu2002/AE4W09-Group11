

%% Rotor Parameters

B = 3; % Number of Blades
R = 143/2; % Insert correct rotor radius
Pitch = 0; % NREL 5 MW has 0 set pitch but this could change

NREL5MW = load("NREL5MW.mat");

%% Operational Specs

lambda = 7.6; % Design tip speed ratio
U0 = 7.5; % U_infinity, average wind speed
rho = 1.225;

%% Airfoil data
Airfoil1.alpha=NREL5MW.Airfoil.Alpha(1);
Airfoil1.Cl=NREL5MW.Airfoil.Cl(1);
Airfoil1.Cd=NREL5MW.Airfoil.Cd(1);

Airfoil2.alpha=NREL5MW.Airfoil.Alpha(2);
Airfoil2.Cl=NREL5MW.Airfoil.Cl(2);
Airfoil2.Cd=NREL5MW.Airfoil.Cd(2);

Airfoil3.alpha=NREL5MW.Airfoil.Alpha(3);
Airfoil3.Cl=NREL5MW.Airfoil.Cl(3);
Airfoil3.Cd=NREL5MW.Airfoil.Cd(3);

Airfoil4.alpha=NREL5MW.Airfoil.Alpha(4);
Airfoil4.Cl=NREL5MW.Airfoil.Cl(4);
Airfoil4.Cd=NREL5MW.Airfoil.Cd(4);

Airfoil5.alpha=NREL5MW.Airfoil.Alpha(5);
Airfoil5.Cl=NREL5MW.Airfoil.Cl(5);
Airfoil5.Cd=NREL5MW.Airfoil.Cd(5);

Airfoil6.alpha=NREL5MW.Airfoil.Alpha(6);
Airfoil6.Cl=NREL5MW.Airfoil.Cl(6);
Airfoil6.Cd=NREL5MW.Airfoil.Cd(6);

Airfoil7.alpha=NREL5MW.Airfoil.Alpha(7);
Airfoil7.Cl=NREL5MW.Airfoil.Cl(7);
Airfoil7.Cd=NREL5MW.Airfoil.Cd(7);

Airfoil8.alpha=NREL5MW.Airfoil.Alpha(8);
Airfoil8.Cl=NREL5MW.Airfoil.Cl(8);
Airfoil8.Cd=NREL5MW.Airfoil.Cd(8);

AirfoilSet = [Airfoil1 Airfoil2 Airfoil3 Airfoil4 Airfoil5 Airfoil6 Airfoil7 Airfoil8];

%% BEM

mu = NREL5MW.Blade.Radius / NREL5MW.Blade.Radius(end);
Twist=NREL5MW.Blade.Twist;
Chord=NREL5MW.Blade.Chord;

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
        if mu(j) < 0.106 % Cylinder
            Airfoil = AirfoilSet(1);
        elseif mu(j) < 0.154 && mu(j) >= 0.106 
            Airfoil = AirfoilSet(2);
        elseif mu(j) < 0.217 && mu(j) >= 0.154
            Airfoil = AirfoilSet(3);
        elseif mu(j) < 0.344 && mu(j) >= 0.217  
            Airfoil = AirfoilSet(4);
        elseif mu(j) < 0.408 && mu(j) >= 0.344   
            Airfoil = AirfoilSet(5);
        elseif mu(j) < 0.535 && mu(j) >= 0.408   
            Airfoil = AirfoilSet(6);
        elseif mu(j) < 0.662 && mu(j) >= 0.535   
            Airfoil = AirfoilSet(7);
        else                                 
            Airfoil = AirfoilSet(8);
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
    FAxial_NREL = FAxialSpan(1,:).*dr; 
    FTang_NREL = FTangSpan(1,:).*dr;

figure
hold on
plot(mu,FAxialSpan,"DisplayName","Axial Force")
plot(mu,FTangSpan,"DisplayName","Tangential Force")
title("Force in each blade section along span, NREL")
xlabel("Spanwise point \mu")
ylabel("Force per unit span [N/m]")
legend
grid on
