function [LPC] = LPCFunc(Ey_kWh,D,Dorig,Pr,Porig)
    %DOESN'T OPTIMIZE FOR D YET
    Dscale = Dorig*((Pr/Porig)^(1/2));              %Is this correct?
    LPC = (0.7+(0.3*((D/Dscale)^2.6)))/Ey_kWh;      %Is this correct?
end 