function [r_list] = gbx_design(r_total, m, types)
    % function that calculates the ratio per stage
    plan_limit = 4;%?????
    para_limit = 7;%?????

    r_approx = nthroot(r_total,m);
    r_list = ones(1,m);

    m_plan = 0;
    for i = 1:m
        if (types(1,i) == true) & (r_approx < plan_limit)
            r_list(1,i) = floor(r_approx);
            m_plan = m_plan + 1;
        elseif (types(1,i) == true)
            r_list(1,i) = plan_limit;
            m_plan = m_plan + 1;
        end
    end

    m_para = m - m_plan;

    r_plans = prod(r_list);
    r_approx = nthroot((r_total/r_plans),m_para);

    for i = 1:m
        if (types(1,i) == false) & (r_approx < para_limit)
            r_list(1,i) = floor(r_approx*10)/10;
        elseif (types(1,i) == false)
            r_list(1,i) = para_limit;
        end
    end

    r_total_res = prod(r_list);
    err = 100*abs(r_total - r_total_res)/r_total;

    disp('The gear ratio was reached with an % error of:')
    disp(err)
end