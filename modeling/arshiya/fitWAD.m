function [results_WAD] = fitWAD(param_struct, data)

numSubj = length(data);

lik = @(x,d) getLogLik(x, d);

logpost_WAD = zeros(numSubj, 1);
loglik_WAD = zeros(numSubj, 1);
best_fit_params_WAD = zeros(numSubj, length(param_struct));
BICs_WAD = zeros(numSubj, 1);
AICs_WAD = zeros(numSubj, 1);

numParams = length(struct);

parfor s = 1:numSubj
    disp(['Fitting WAD for subject ', num2str(s)]);

    WAD_post = @(x) -(lik(x,data(s)) + getPriorSum(x, param_struct));
    [x,logpost] = ga(WAD_post, length(param_struct),[],[],[],[],vertcat(param_struct.lb), ...
        vertcat(param_struct.ub),[],find(vertcat(param_struct.int)));
    best_fit_params_WAD(s,:) = x;
    loglik_WAD(s) = lik(x, data(s));
    logpost_WAD(s) = -logpost;
    BICs_WAD(s) = numParams*log(data(s).N) - 2*loglik_WAD(s);
    AICs_WAD(s) = numParams*2 - 2*loglik_WAD(s);

    disp(['Completed optimization for WAD for subject ', num2str(s)]);
end

results_template.K = numParams;
results_template.S = numSubj;
results_template.param = param_struct;

results_WAD = results_template;
results_WAD.logpost = logpost_WAD';
results_WAD.loglik = loglik_WAD;
results_WAD.x = best_fit_params_WAD;
results_WAD.bic = BICs_WAD;
results_WAD.aic = AICs_WAD;
results_WAD.likfun = lik;

end