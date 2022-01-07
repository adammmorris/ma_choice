function priorsum = getPriorSum(params, param_struct)
priorsum = 0;
for i = 1:length(param_struct)
    priorsum = priorsum + param(i).logpdf(params(i));
end