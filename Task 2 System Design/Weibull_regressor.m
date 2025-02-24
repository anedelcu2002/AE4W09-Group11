function [f_curve] = Weibull_regressor(U_array) 
    % first, transforms the wind speed time-series U_array in m/s into a histogram, then fits a Weibull distribution onto that histogram

    U_array=transpose(U_array);
    pd=fitdist(U_array,'wbl');
    parameters=paramci(pd);
    a=(parameters(1,1)+parameters(2,1))/2;
    k=(parameters(1,2)+parameters(2,2))/2;


    f_curve = [];

    for i = 1:ceil(max(U_array))
        f_curve(i) = (k/a)*((i/a)^(k-1))*exp(-(i/a)^k);
    end

end