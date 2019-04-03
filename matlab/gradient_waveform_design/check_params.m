function Algo_param = check_params(Algo_param)
% function Algo_param = check_params(Algo_param)
% check parameters and set default ones if needed.

if ~isfield(Algo_param,'show_progression')
    Algo_param.show_progression=0;
end

if ~isfield(Algo_param,'display_results')
    Algo_param.display_results=0;
end

