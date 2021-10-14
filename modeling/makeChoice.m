function [choices] = makeChoice(parameters, options, avail_atts)

% parameters: [invTemp weights]
% options: numAtt x numOptionsPerChoice x numChoices
% avail_atts: numChoices x numAtt (logicals)

options_dim = size(options);

if nargin < 3 || isempty(avail_atts); avail_atts = true(options_dim(3), options_dim(1)); end
choices = zeros(options_dim(3), 1);
inv_temp = parameters(1);
weights = parameters(2:end);

% compute utilities & make choices
for i = 1:options_dim(3)
    utilities = weights(avail_atts(i,:)) * options(avail_atts(i,:),:,i);
    probs = exp(inv_temp * utilities - logsumexp(inv_temp * utilities, 2));
    choices(i) = fastrandsample(probs);
end

end