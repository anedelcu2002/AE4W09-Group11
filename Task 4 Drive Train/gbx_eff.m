function [eff_gbx] = gbx_eff(r_list, g_types)
    eff_list1 = [];
    eff_list2 = [];
    for i = 1:length(r_list)
        if g_types(i) == true
            eff0 = 0.97;
        elseif g_types(i) == false
            eff0 = 0.98;
        end

        eff_list1(i) = eff0 * exp(-0.001 * (r_list(i)-1));
        eff_list2(i) = eff0;
    end

    %eff_gbx = prod(eff_list1);
    eff_gbx = prod(eff_list2);
end