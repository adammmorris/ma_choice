function [loglik] = getLogLik(parameters, data, sign_optionvals)

% parameters: invTemp weights
% data: structure with .choices and .options

options = signOptions(data.options, sign_optionvals, true);
choices = data.choices;
options_dim = size(options);
numAtts = options_dim(1);
numChoices = options_dim(3);
if isfield(data, 'avail_atts')
    avail_atts = logical(data.avail_atts);
else
    avail_atts = true(numChoices, numAtts);
end
loglik = 0;
inv_temp = parameters(1);
weights = parameters(2:end);

for i = 1:numChoices
    utilities = weights(avail_atts(i,:)) * options(avail_atts(i,:),:,i);
    loglik = loglik + inv_temp * utilities(choices(i)) - logsumexp(inv_temp * utilities, 2);
end

end