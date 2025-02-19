%Original Turbine
Porig = 5*10^6;
Dorig = 125;%??????

%Constants & Pre-decided values
k = 2;
a_scale = 8; 
Pr = 7.5*10^6;
rho = 1.225;

%Approximate values
eff = 0.9;
Cp = 0.5;
%%
%Optimize R
Uci = 5;
Uco = 30;
R_range = linspace(0,200,200);

Ey_list = [];
LPC_list = [];

for i = 1:length(R_range)
    Ey_list(i) = EnergyFunc(R_range(i),Uci,Uco,k,a_scale,Pr,rho,eff,Cp);
    LPC_list(i) = LPCFunc(Ey_list(i),2*R_range(i),Dorig,Pr,Porig);
end

figure
plot(R_range,Ey_list)
grid on

figure
plot(R_range,LPC_list)
grid on
%%
R = 80; %Choice based on above figures
%%
%Optimize Uci,Uco
Uci_range = linspace(0,100,100);
Uco_range = linspace(0,100,100);

Ey_matrix = [];
LPC_matrix = [];

for i = 1:length(Uci_range)
    for j = 1:length(Uco_range)
        if Uco_range(j) > Uci_range(i)
            Ey_matrix(j,i) = EnergyFunc(R,Uci_range(i),Uco_range(j),k,a_scale,Pr,rho,eff,Cp);
        else
            Ey_matrix(j,i) = 0;
        end
        LPC_matrix(j,i) = LPCFunc(Ey_matrix(j,i),2*R,Dorig,Pr,Porig);
    end
end

figure
set(gcf, 'Position',[10 10 800 500])
surf(Uci_range,Uco_range, Ey_matrix,'FaceColor','white','FaceAlpha',0.9); hold on;
view(3); % Set 3D view
xlabel('Uci');ylabel('Uco');zlabel('Ey');grid on;

figure
set(gcf, 'Position',[10 10 800 500])
surf(Uci_range,Uco_range, LPC_matrix,'FaceColor','white','FaceAlpha',0.9); hold on;
view(3); % Set 3D view
xlabel('Uci');ylabel('Uco');zlabel('LPC');grid on;
%%
%Choice based on above figures
Uci = 5;
Uco = 25;
%%
Ur = (Pr/(0.5*eff*rho*Cp*pi*R*R))^(1/3)

z0 = 0.03;
alpha = 0.143;

href = 60; %should be between 60 and 100?
Uref = 13.1;


%Hub height
h_list = linspace(0, 100, 100);
Uh_list = [];

for n = 1:length(h_list)
    if h_list(n) < 60;
        Uh_list(n) = Uref*((log(h_list(n)/z0)) / (log(href/z0)));
    else
        Uh_list(n) = Uref*((h_list(n)/href)^alpha);
    end
end

for n = 1:length(h_list)
    if Uh_list(n) > Ur
        h = h_list(n)
        break
    end
end


figure
hold on
plot(Uh_list,h_list)
yline(h)
xline(Ur)
hold off
grid on
%%
lambda = 7; %To be chosen???

Omega_r = (lambda*Ur)/R %rad/s
V_tip = Omega_r*R  %m/s
%%
%Aerodynamic Torque --> Torque on LSS
Qlss = Pr/Omega_r