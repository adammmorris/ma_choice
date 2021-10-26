function [results_WAD, results_WP, results_EW, results_TAL] = fitModels(param, data, nstarts, numAtts)

lbs = [param(1).lb repelem(-1, numAtts)];
ubs = [param(1).ub repelem(1, numAtts)];
numSubj = length(data);
hessians = cell(numSubj, 1);
numParams = length(param);

results_WAD = mfit_optimize_parallel(@getLogLik_WAD,param,data,nstarts);
results_WP = mfit_optimize_parallel(@getLogLik_WP,param,data,nstarts);

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

% logpost_LEX = zeros(numSubj, 1);
% loglik_LEX = zeros(numSubj, 1);
% best_fit_params_LEX = zeros(numSubj, numParams);
% BICs_LEX = zeros(numSubj, 1);
% AICs_LEX = zeros(numSubj, 1);

parfor s = 1:numSubj
    % EW
    EW_post = @(x) -(getLogLik_WAD(x,data(s)) + param(1).logpdf(x(1)) + log((1/3) ^ numAtts));

    [x,logpost] = ga(EW_post, numParams,[],[],[],[],lbs,ubs,[],2:numParams);
    best_fit_params_EW(s,:) = x;
    loglik_EW(s) = getLogLik_WAD(x, data(s));
    logpost_EW(s) = -logpost;
    BICs_EW(s) = numParams*log(data(s).N) - 2*loglik_EW(s);
    AICs_EW(s) = numParams*2 - 2*loglik_EW(s);

    % TAL
    TAL_post = @(x) -(getLogLik_WP(x,data(s)) + param(1).logpdf(x(1)) + log((1/3) ^ numAtts));

    [x,logpost] = ga(TAL_post, numParams,[],[],[],[],lbs,ubs,[],2:numParams);
    best_fit_params_TAL(s,:) = x;
    loglik_TAL(s) = getLogLik_WP(x, data(s));
    logpost_TAL(s) = -logpost;
    BICs_TAL(s) = numParams*log(data(s).N) - 2*loglik_TAL(s);
    AICs_TAL(s) = numParams*2 - 2*loglik_TAL(s);

    % LEX
%     LEX_post = @(x) -(getLogLik_WAD(x,data(s)) + param(1).logpdf(x(1)) + ...
%         log(1 / (numAtts * 2)));
%     
%     [x,logpost] = ga(LEX_post, numParams,...
%         [],[],[],[],...
%         lbs,ubs,@my_nonlcon,2:numParams);
%     best_fit_params_LEX(s,:) = x;
%     loglik_LEX(s) = getLogLik_WP(x, data(s));
%     logpost_LEX(s) = -logpost;
%     BICs_LEX(s) = numParams*log(data(s).N) - 2*loglik_LEX(s);
%     AICs_LEX(s) = numParams*2 - 2*loglik_LEX(s);
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

results_TAL = results_WP;
results_TAL.K = numParams;
results_TAL.S = numSubj;
results_TAL.H = hessians;
results_TAL.logpost = logpost_TAL';
results_TAL.loglik = loglik_TAL;
results_TAL.x = best_fit_params_TAL;
results_TAL.bic = BICs_TAL;
results_TAL.aic = AICs_TAL;

% results_LEX = results_WAD;
% results_LEX.K = numParams;
% results_LEX.S = numSubj;
% results_LEX.H = hessians;
% results_LEX.logpost = logpost_LEX';
% results_LEX.loglik = loglik_LEX;
% results_LEX.x = best_fit_params_LEX;
% results_LEX.bic = BICs_LEX;
% results_LEX.aic = AICs_LEX;

end