function [eff_gbx] = gbx_eff(r_list, g_types)
    eff_list = [];
    for i = 1:length(r_list)
        if g_types(i) == true
            eff0 = 0.99;
        elseif g_types(i) == false
            eff0 = 0.98;
        end
        eff_list(i) = eff0;
    end

    eff_gbx = prod(eff_list);
end