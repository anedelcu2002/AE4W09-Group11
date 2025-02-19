function [f_curve] = WeibullFunc(k,a_scale,U_range) 
%THIS ONE IS STILL BASED ON k AND a

%Weibelcurve
f_curve = [];

for i = linspace(1,length(U_range),length(U_range))
    f_curve(i) = (k/a_scale)*((U_range(i)/a_scale)^(k-1))*exp(-(U_range(i)/a_scale)^k);
end

end