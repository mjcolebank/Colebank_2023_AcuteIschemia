% Function to reduce active force as a function of time

function percent_reduce = induce_MI(pars,t)
NC_stopMI = pars.MI(1);
NC_induce_MI =  pars.MI(2);
t_end = pars.T*NC_stopMI;
t_MI = pars.T*NC_induce_MI;
act_max = 1.0;
act_min = pars.MI(3);
if isscalar(t)
    if t<t_MI
        percent_reduce = 1.0;
    elseif t<t_end
        percent_reduce = 1.0+(t-t_MI).*(act_min-act_max)./(t_end-t_MI);
%             percent_reduce
    else
        percent_reduce = act_min;
    end
else
    percent_reduce = 1.0+max(0,t-t_MI).*(act_min-act_max)./(t_end-t_MI);
    ids = t>=t_end;
    percent_reduce(ids) = act_min;
end

end