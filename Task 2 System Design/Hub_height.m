function [h_hub] = Hub_height(D,D_orig,h_hub_orig)
    % scale hub height linearly with diameter based on original turbine

    h_hub=h_hub_orig*(D/D_orig);     
end 