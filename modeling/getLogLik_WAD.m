function [loglik] = getLogLik_WAD(parameters, data)

% parameters: invTemp weights
% data: structure with .choices and .options

options = data.options;
choices = data.choices;
options_dim = size(options);
loglik = 0;
inv_temp = parameters(1);
weights = parameters(2:end);

% compute utilities & make choices
for i = 1:options_dim(3)
    utilities = weights * options(:,:,i);
    loglik = loglik + inv_temp * utilities(choices(i)) - logsumexp(inv_temp * utilities, 2);
end

end