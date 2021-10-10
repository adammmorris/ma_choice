function [options_signed] = signOptions(options)
% options should be numAtts x 2 x numChoices x numAgents
[numAtts, ~, numChoices, numAgents] = size(options);
options_signed = zeros(numAtts, 2, numChoices, numAgents);
for i = 1:numAgents
    for j = 1:numChoices
        options_signed(:,1,j,i) = sign(options(:,1,j,i) - options(:,2,j,i));
    end
end
end