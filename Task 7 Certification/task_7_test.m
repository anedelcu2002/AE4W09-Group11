%% to-do:
% calculate e-modulus based on literature
% calculate damage equivalent load for fatigue
% switch stress calculation to geometric vector addition, variable max stress point, distance to neutral axis is distance on circle!!
% justify use of superimposed loads for fatigue, reflect on how it is an improper assessment of amplitudes
% check whether to use moments or stresses for fatigue
% add safety factors for deflection

function cumulativeDamageStemPlot(ni,Nfi) % stolen from https://nl.mathworks.com/help/signal/ug/practical-introduction-to-fatigue-analysis-using-rainflow-counting.html
figure
L = length(ni);
damage = sum(ni./Nfi);
stem(0,NaN,"Color",[0 1 0])
title("Cumulative Damage from Palmgren-Miner Rule")
xlabel("Cycle $i$","Interpreter","latex")
ylabel("Cum. damage $D_{i} = \sum_{j=1}^{i}n_{j}/N_{f,j}$","Interpreter","latex")
set(gca,"XLim",[0 L],"YLim",[0 damage])
grid("on")
iter = unique([1:round(L/100):L,L]);
hold(gca,"on")
for i = iter
    cdi = sum(ni(1:i)./Nfi(1:i)); % cumulative damge upto cycle i
    plt = stem(i,cdi,"filled");
    setStemColor(plt,cdi,0.95)
end
end

function setStemColor(hplt,cumulativeDamage,gamma)
c = lines(5);
c = c([2,3,5],:);
if (cumulativeDamage > 1)
    color = c(1,:);
else
    if (cumulativeDamage > gamma)
        c1 = c(1,:);
        c2 = c(2,:);
    else
        c1 = c(3,:);
        c2 = c(2,:);
    end
    color = zeros(1,3);    
    for i = 1:3
        color(i) = c1(i)+(c2(i)-c1(i))*cumulativeDamage;
    end
end
hplt.Color = color;
end

%% initialize data
response_data=load('nrel_cert.mat');
response_data.Legend;

turbine_data=load('NREL5MW.mat');

E=4.7*10^9;
EI_edge=turbine_data.Blade.EIedge(1);
EI_flap=turbine_data.Blade.EIflap(1);
thickness=turbine_data.Blade.Thickness(1)/sqrt(2); %base root is a cylinder, assuming spar structure is square enclosed in circle

%% perform post-processing
analysis=0; % 1 for extreme load, 0 for fatigue

if analysis==1

    figure;
    plot(response_data.Time, response_data.OoPDefl1);
    hold on;
    plot(response_data.Time, response_data.OoPDefl2);
    hold on;
    plot(response_data.Time, response_data.OoPDefl3);
    legend('blade 1', 'blade 2', 'blade 3');
    title('out of plane tip deflection')
    ylabel('deflection (m)')
    xlabel('time (s)')

    figure;
    plot(response_data.Time, response_data.RootMEdg1);
    hold on;
    plot(response_data.Time, response_data.RootMEdg2);
    hold on;
    plot(response_data.Time, response_data.RootMEdg3);
    legend('blade 1', 'blade 2', 'blade 3');
    title('blade root moment edge direction')
    ylabel('moment (kN m)')
    xlabel('time (s)')

    figure;
    plot(response_data.Time, response_data.RootMFlp1);
    hold on;
    plot(response_data.Time, response_data.RootMFlp2);
    hold on;
    plot(response_data.Time, response_data.RootMFlp3);
    legend('blade 1', 'blade 2', 'blade 3');
    title('blade root moment flap direction')
    ylabel('moment (kN m)')
    xlabel('time (s)')

    maximum_deflection=max([max(response_data.OoPDefl1), max(response_data.OoPDefl2), max(response_data.OoPDefl3)]); %maximum out of plane tip deflection in meters

    root_stress1=(response_data.RootMFlp1*1000*thickness/2/EI_flap+response_data.RootMEdg1*1000*thickness/2/EI_edge)*E; %blade 1 root stress timeseries
    root_stress2=(response_data.RootMFlp2*1000*thickness/2/EI_flap+response_data.RootMEdg2*1000*thickness/2/EI_edge)*E; 
    root_stress3=(response_data.RootMFlp3*1000*thickness/2/EI_flap+response_data.RootMEdg3*1000*thickness/2/EI_edge)*E;

    gam_m=1.35; % safety factor for material properties
    gam_f=1.3; % safety factor for loads
    gam_n=1; % safety factor for severity of effect
    maximum_stress=gam_m*gam_f*gam_n*max([max(root_stress1), max(root_stress2), max(root_stress3)]); %maximum blade root stress
    maximum_stress_mpa=maximum_stress/10^6

elseif analysis==0

    figure;
    plot(response_data.Time, response_data.OoPDefl1);
    hold on;
    plot(response_data.Time, response_data.OoPDefl2);
    hold on;
    plot(response_data.Time, response_data.OoPDefl3);
    legend('blade 1', 'blade 2', 'blade 3');

    figure;
    plot(response_data.Time, response_data.RootMEdg1);
    hold on;
    plot(response_data.Time, response_data.RootMEdg2);
    hold on;
    plot(response_data.Time, response_data.RootMEdg3);
    legend('blade 1', 'blade 2', 'blade 3');
    title('blade root moment edge direction')
    ylabel('moment (kN m)')
    xlabel('time (s)')

    figure;
    plot(response_data.Time, response_data.RootMFlp1);
    hold on;
    plot(response_data.Time, response_data.RootMFlp2);
    hold on;
    plot(response_data.Time, response_data.RootMFlp3);
    legend('blade 1', 'blade 2', 'blade 3');
    title('blade root moment flap direction')
    ylabel('moment (kN m)')
    xlabel('time (s)')

    root_stress1=(response_data.RootMFlp1*1000*thickness/2/EI_flap+response_data.RootMEdg1*1000*thickness/2/EI_edge)*E; %blade 1 root stress timeseries
    root_stress2=(response_data.RootMFlp2*1000*thickness/2/EI_flap+response_data.RootMEdg2*1000*thickness/2/EI_edge)*E; 
    root_stress3=(response_data.RootMFlp3*1000*thickness/2/EI_flap+response_data.RootMEdg3*1000*thickness/2/EI_edge)*E;
    
    figure;
    plot(response_data.Time, root_stress1);
    hold on;
    plot(response_data.Time, root_stress2);
    hold on;
    plot(response_data.Time, root_stress3);
    legend('blade 1', 'blade 2', 'blade 3');
    title('blade root stress amplitude')
    ylabel('stress (N/m2)')
    xlabel('time (s)')

    %% rainflow counting (from brightspace)
    gam_m=1; % safety factor for material properties
    gam_f=1.2; % safety factor for loads
    gam_n=1.15; % safety factor for severity of effect
    S_s1 = gam_m*gam_f*gam_n*root_stress1;
    S_s2 = gam_m*gam_f*gam_n*root_stress2;
    S_s3 = gam_m*gam_f*gam_n*root_stress3;

    % Perform rainflow count, blade 1
    [c1,hist1,edges1,rmm1,idx1] = rainflow(S_s1);

    % Plot histogram
    figure;
    histogram('BinEdges',edges1','BinCounts',sum(hist1,2))
    title('Blade 1 rainflow counting')
    xlabel('Stress Range')
    ylabel('Cycle Counts')

    % Set S-N curve and use Minor's rule to determine damage
    %   by summing damages from each individual cycle or half cycle
    m=9; % slope S-N curve N*S^m=K (with S stress range)
    UCS=600*10^6; % ultimate compression strength (Nm2);
    K=(2*UCS)^m; % S-N in stress ranges, so factor 2 for amplitudes

    range= c1(:,2);
    count= c1(:,1);
    D1=1/K*sum(count.*(range.^m));

    cumulativeDamageStemPlot(count,K./(range.^m));

    % Perform rainflow count, blade 2
    [c2,hist2,edges2,rmm2,idx2] = rainflow(S_s2);

    figure;
    histogram('BinEdges',edges2','BinCounts',sum(hist2,2))
    title('Blade 2 rainflow counting')
    xlabel('Stress Range')
    ylabel('Cycle Counts')

    range= c2(:,2);
    count= c2(:,1);
    D2=1/K*sum(count.*(range.^m));

    cumulativeDamageStemPlot(count,K./(range.^m));

    % Perform rainflow count, blade 3
    [c3,hist3,edges3,rmm3,idx3] = rainflow(S_s3);

    figure;
    histogram('BinEdges',edges3','BinCounts',sum(hist3,2))
    title('Blade 3 rainflow counting')
    xlabel('Stress Range')
    ylabel('Cycle Counts')

    range= c3(:,2);
    count= c3(:,1);
    D3=1/K*sum(count.*(range.^m));

    cumulativeDamageStemPlot(count,K./(range.^m));
end


%M=response_data.RootMFlp1; % blade root flap moment (kNm)
% Remove 60 s start-up transient
%Ind=find(response_data.Time<=60);
%M(Ind)=[];
% Set geometrical properties of the cross section and determine stress 
%I=1.37; % blade root inertia (m^4)
%y=1.5; % distance to neutral line (m)
%S=M*y/I;
% Set safety factors and apply these to stress
%gam_m=1; % safety factor for material properties
%gam_f=1; % safety factor for loads
%gam_n=1; % safety factor for severity of effect
%S_s = gam_m*gam_f*gam_n*S;

% Perform rainflow count
%[c,hist,edges,rmm,idx] = rainflow(S_s);

% Plot histogram
%histogram('BinEdges',edges','BinCounts',sum(hist,2))
%xlabel('Stress Range')
%ylabel('Cycle Counts')

% Set S-N curve and use Minor's rule to determine damage
%   by summing damages from each individual cycle or half cycle
%m=9; % slope S-N curve N*S^m=K (with S stress range)
%UCS=6e5; % ultimate compression strength (kNm);
%K=(2*UCS)^m; % S-N in stress ranges, so factor 2 for amplitudes
%range= c(:,2);
%count= c(:,1);
%D=1/K*sum(count.*(range.^m));