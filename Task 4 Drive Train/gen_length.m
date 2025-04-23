function [l_s] = gen_length(V_gen, r_s)
    % function that calculates the generator length
        l_s = V_gen/(pi*r_s*r_s); 
end