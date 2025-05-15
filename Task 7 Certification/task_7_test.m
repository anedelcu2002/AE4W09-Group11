%% to-do:
% calculate e-modulus based on literature - alex
% determine wind speeds for each simulation - alex
% justify use of superimposed loads for fatigue, reflect on how it is an improper assessment of amplitudes - alex
% document all work done so far - alex

% calculate damage equivalent load for fatigue - jesse
% check simulation duration for fatigue - jesse


% plot some nice spectra - victor


% calculate thickness factor and recheck design with task 3, task 5 - victor

function cumulativeDamageStemPlot(ni,Nfi, blade_nr) % stolen from https://nl.mathworks.com/help/signal/ug/practical-introduction-to-fatigue-analysis-using-rainflow-counting.html
    figure
    L = length(ni);
    damage = sum(ni./Nfi);
    stem(0,NaN,"Color",[0 1 0])
    title(["Cumulative Damage from Palmgren-Miner Rule, blade ", num2str(blade_nr)])
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

function stress=calculate_stress(moment_1, moment_2, thickness, EI_1, EI_2, E)
    moment_sum=(moment_1.^2.+moment_2.^2).^0.5;
    distance_1=thickness.*(moment_1./moment_sum);
    distance_2=thickness.*(moment_2./moment_sum);
    stress=moment_1.*distance_1./EI_1.*E+moment_2.*distance_2./EI_2.*E;
end

function D = rainflow_counting(stress_timeseries, blade_nr, plot_b) % stolen from brightspace
    gam_m=1; % safety factor for material properties
    gam_f=1.2; % safety factor for loads
    gam_n=1.15; % safety factor for severity of effect

    S_s= gam_m.*gam_f.*gam_n.*stress_timeseries;

    % Perform rainflow count, blade 1
    [c,hist,edges,rmm,idx] = rainflow(S_s);

    % Plot histogram
    if plot_b
        figure;
        histogram('BinEdges',edges','BinCounts',sum(hist,2))
        title(['Rainflow counting, blade ', num2str(blade_nr)])
        xlabel('Stress Range')
        ylabel('Cycle Counts')
    end

    % Set S-N curve and use Minor's rule to determine damage
    %   by summing damages from each individual cycle or half cycle
    m=9; % slope S-N curve N*S^m=K (with S stress range)
    UCS=600*10^6; % ultimate compression strength (Nm2);
    K=(2*UCS)^m; % S-N in stress ranges, so factor 2 for amplitudes

    range= c(:,2);
    count= c(:,1);
    D=1/K*sum(count.*(range.^m));
    if plot_b
        cumulativeDamageStemPlot(count,K./(range.^m), blade_nr);
    end
end


%% initialize data
response_data=load('nrel_cert.mat');
response_data.Legend;

turbine_data=load('NREL5MW.mat');

E=4.7*10^9;
EI_edge=turbine_data.Blade.EIedge(1);
EI_flap=turbine_data.Blade.EIflap(1);
thickness=turbine_data.Blade.Thickness(1); %base root is a load-bearing cylinder

%% plot simulation outputs
plot_b=true;

if plot_b
    figure;
    plot(response_data.Time, response_data.OoPDefl1);
    hold on;
    plot(response_data.Time, response_data.OoPDefl2);
    hold on;
    plot(response_data.Time, response_data.OoPDefl3);
    legend('blade 1', 'blade 2', 'blade 3');
    title('Out of plane tip deflection')
    ylabel('deflection (m)')
    xlabel('time (s)')

    figure;
    plot(response_data.Time, response_data.RootMEdg1);
    hold on;
    plot(response_data.Time, response_data.RootMEdg2);
    hold on;
    plot(response_data.Time, response_data.RootMEdg3);
    legend('blade 1', 'blade 2', 'blade 3');
    title('Blade root moment, edge direction')
    ylabel('moment (kN m)')
    xlabel('time (s)')

    figure;
    plot(response_data.Time, response_data.RootMFlp1);
    hold on;
    plot(response_data.Time, response_data.RootMFlp2);
    hold on;
    plot(response_data.Time, response_data.RootMFlp3);
    legend('blade 1', 'blade 2', 'blade 3');
    title('Blade root moment, flap direction')
    ylabel('moment (kN m)')
    xlabel('time (s)')
end

%% perform post-processing
analysis=0; % 1 for extreme load, 0 for fatigue

if analysis==1
    maximum_deflection=max([max(response_data.OoPDefl1), max(response_data.OoPDefl2), max(response_data.OoPDefl3)]); %maximum out of plane tip deflection in meters
    sf_m_defl=1.1; % safety factor for material properties
    sf_f_defl=1.35; % safety factor for loads
    sf_n_defl=1; % safety factor for severity of effect
    maximum_deflection=maximum_deflection*sf_m_defl*sf_f_defl*sf_n_defl; % maximum deflection in meters
    disp(['Maximum tip deflection in out-of-plane direction is ', num2str(maximum_deflection), ' meters'])

    root_stress1=calculate_stress(response_data.RootMFlp1*1000, response_data.RootMEdg1*1000, thickness/2, EI_flap, EI_edge, E);%blade 1 root stress timeseries
    root_stress2=calculate_stress(response_data.RootMFlp2*1000, response_data.RootMEdg2*1000, thickness/2, EI_flap, EI_edge, E);
    root_stress3=calculate_stress(response_data.RootMFlp3*1000, response_data.RootMEdg3*1000, thickness/2, EI_flap, EI_edge, E);

    if plot_b
        figure;
        plot(response_data.Time, root_stress1);
        hold on;
        plot(response_data.Time, root_stress2);
        hold on;
        plot(response_data.Time, root_stress3);
        legend('blade 1', 'blade 2', 'blade 3');
        title('Blade root stress amplitude')
        ylabel('stress (N/m2)')
        xlabel('time (s)')
    end

    gam_m=1.35; % safety factor for material properties
    gam_f=1.3; % safety factor for loads
    gam_n=1; % safety factor for severity of effect

    maximum_stress=gam_m*gam_f*gam_n*max([max(root_stress1), max(root_stress2), max(root_stress3)]); %maximum blade root stress after safety factors
    disp(['Maximum blade root stress is ', num2str(maximum_stress/10^6), ' MPa'])

elseif analysis==0
    root_stress1=calculate_stress(response_data.RootMFlp1*1000, response_data.RootMEdg1*1000, thickness/2, EI_flap, EI_edge, E);%blade 1 root stress timeseries
    root_stress2=calculate_stress(response_data.RootMFlp2*1000, response_data.RootMEdg2*1000, thickness/2, EI_flap, EI_edge, E);
    root_stress3=calculate_stress(response_data.RootMFlp3*1000, response_data.RootMEdg3*1000, thickness/2, EI_flap, EI_edge, E);

    if plot_b
        figure;
        plot(response_data.Time, root_stress1);
        hold on;
        plot(response_data.Time, root_stress2);
        hold on;
        plot(response_data.Time, root_stress3);
        legend('blade 1', 'blade 2', 'blade 3');
        title('Blade root stress amplitude')
        ylabel('stress (N/m2)')
        xlabel('time (s)')
    end

    D1 = rainflow_counting(root_stress1, 1, plot_b);
    D2 = rainflow_counting(root_stress2, 2, plot_b);
    D3 = rainflow_counting(root_stress3, 3, plot_b);
    disp(['Cumulative fatigue damage on blade 1 is ', num2str(D1)])
    disp(['Cumulative fatigue damage on blade 2 is ', num2str(D2)])
    disp(['Cumulative fatigue damage on blade 3 is ', num2str(D3)])

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