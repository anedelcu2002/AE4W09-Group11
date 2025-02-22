function [f_curve] = Weibull_regressor(U_array) 
    % first, transforms the wind speed time-series U_array in m/s into a histogram, then fits a Weibull distribution onto that histogram

    U_sorted=histogram(U_array, max(U_array)); %number of bins is equal to the maximum wind speed, so that the bins are split by 1 m/s

    [k, a]=wblfit(U_sorted); %shape and scale parameter

    f_curve = [];

    for i = linspace(1,length(U_range),max(U_array)-min(U_array))
        f_curve(i) = (k/a)*((U_range(i)/a)^(k-1))*exp(-(U_range(i)/a)^k);
    end

end