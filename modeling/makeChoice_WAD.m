function [choices] = makeChoice_WAD(parameters, options)

% parameters: [invTemp weights]
% options: numAtt x numOptionsPerChoice x numChoices

options_dim = size(options);
choices = zeros(options_dim(3), 1);
inv_temp = parameters(1);
weights = parameters(2:end);

% compute utilities & make choices
for i = 1:options_dim(3)
    utilities = weights * options(:,:,i);
    probs = exp(inv_temp * utilities - logsumexp(inv_temp * utilities, 2));
    choices(i) = fastrandsample(probs);
end

end