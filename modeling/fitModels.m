function [results_WAD, results_WP, results_EW, results_TAL] = fitModels(param_struct, param_struct_signedweights, data, nstarts)

numSubj = length(data);
hessians = cell(numSubj, 1);

lik = @(x,d) getLogLik(x, d, false);
lik_signedopts = @(x,d) getLogLik(x, d, true);

results_WAD = mfit_optimize_parallel(lik,param_struct,data,nstarts);
results_WP = mfit_optimize_parallel(lik_signedopts,param_struct,data,nstarts);

logpost_EW = zeros(numSubj, 1);
loglik_EW = zeros(numSubj, 1);
best_fit_params_EW = zeros(numSubj, numParams);
BICs_EW = zeros(numSubj, 1);
AICs_EW = zeros(numSubj, 1);

logpost_TAL = zeros(numSubj, 1);
loglik_TAL = zeros(numSubj, 1);
best_fit_params_TAL = zeros(numSubj, numParams);
BICs_TAL = zeros(numSubj, 1);
AICs_TAL = zeros(numSubj, 1);

parfor s = 1:numSubj
    % EW
    EW_post = @(x) -(lik(x,data(s)) + getPriorSum(x, param_struct_signedweights));
    [x,logpost] = ga(EW_post, length(param_struct_signedweights),[],[],[],[],vertcat(param_struct_signedweights.lb), ...
        vertcat(param_struct_signedweights.ub),[],2:length(param_struct_signedweights));
    best_fit_params_EW(s,:) = x;
    loglik_EW(s) = lik(x, data(s));
    logpost_EW(s) = -logpost;
    BICs_EW(s) = numParams*log(data(s).N) - 2*loglik_EW(s);
    AICs_EW(s) = numParams*2 - 2*loglik_EW(s);

    % TAL
    TAL_post = @(x) -(lik_signedopts(x,data(s)) + getPriorSum(x, param_struct_signedweights));
    [x,logpost] = ga(TAL_post, length(param_struct_signedweights),[],[],[],[],vertcat(param_struct_signedweights.lb), ...
        vertcat(param_struct_signedweights.ub),[],2:length(param_struct_signedweights));
    best_fit_params_TAL(s,:) = x;
    loglik_TAL(s) = lik_signedopts(x, data(s));
    logpost_TAL(s) = -logpost;
    BICs_TAL(s) = numParams*log(data(s).N) - 2*loglik_TAL(s);
    AICs_TAL(s) = numParams*2 - 2*loglik_TAL(s);
end

results_EW = results_WAD;
results_EW.K = numParams;
results_EW.S = numSubj;
results_EW.H = hessians;
results_EW.logpost = logpost_EW';
results_EW.loglik = loglik_EW;
results_EW.x = best_fit_params_EW;
results_EW.bic = BICs_EW;
results_EW.aic = AICs_EW;
results_EW.likfun = lik;

results_TAL = results_WP;
results_TAL.K = numParams;
results_TAL.S = numSubj;
results_TAL.H = hessians;
results_TAL.logpost = logpost_TAL';
results_TAL.loglik = loglik_TAL;
results_TAL.x = best_fit_params_TAL;
results_TAL.bic = BICs_TAL;
results_TAL.aic = AICs_TAL;
results_TAL.likfun = lik_signedopts;

end