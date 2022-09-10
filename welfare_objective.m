function [outvalue]=welfare_objective(x_opt,x_opt_name)
% function [outvalue]=welfare_objective(x_opt,x_opt_name)

global oo_ options_ M_

%% set parameter for use in Dynare
for ii=1:size(x_opt_name,1)
    set_param_value(x_opt_name{ii,1},x_opt(ii));
end

if any(x_opt<cell2mat(x_opt_name(:,2))) || any(x_opt>cell2mat(x_opt_name(:,3))) %make sure parameters are inside their bounds
    outvalue=10e6+sum((x_opt).^2); %penalty function
    return
end

if isempty(options_.qz_criterium)
    options_.qz_criterium = 1+1e-6;
end

var_list_ = char('welf');
[oo_.dr,info,M_,oo_] = resol(0,M_,options_,oo_); %get decision rules
if info(1) %filter out error code
    outvalue=1e5+sum((x_opt).^2);
    return;
end

%% simulate conditional welfare
initial_condition_states = repmat(oo_.dr.ys,1,M_.maximum_lag); %get steady state as initial condition
shock_matrix = zeros(1,M_.exo_nbr); %create shock matrix with number of time periods in columns
y_sim = simult_(M_,options_,initial_condition_states,oo_.dr,shock_matrix,options_.order); %simulate one period to get value

outvalue=-y_sim(strmatch('welf',M_.endo_names,'exact'),2); %extract Welfare