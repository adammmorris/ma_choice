function [results_WADs, results_WPs, results_EWs, results_TALs] = ...
    fitModels_step(param_struct_step, param_struct_signedweights_step, data, which_models)%, fittedpars, truepars, probs_gen)

if nargin < 5, which_models = {'WAD', 'WP', 'EW', 'TAL'}; end

numSubj = length(data);
hessians = cell(numSubj, 1);
numParams = length(param_struct_step);

logpost_WAD = zeros(numSubj, 1);
loglik_WAD = zeros(numSubj, 1);
best_fit_params_WAD = zeros(numSubj, numParams);
BICs_WAD = zeros(numSubj, 1);
AICs_WAD = zeros(numSubj, 1);

logpost_WP = zeros(numSubj, 1);
loglik_WP = zeros(numSubj, 1);
best_fit_params_WP = zeros(numSubj, numParams);
BICs_WP = zeros(numSubj, 1);
AICs_WP = zeros(numSubj, 1);

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
 
lik = @(x,d) getLogLik(x, d, false);
lik_signedopts = @(x,d) getLogLik(x, d, true);

parfor s = 1:numSubj
    if any(strcmp(which_models, 'WAD'))
        % WAD
        disp(['Fitting WAD_step for subject ', num2str(s)]);
        WAD_post = @(x) -(lik(x,data(s)) + getPriorSum(x, param_struct_step));

%         chance = fittedpars; chance(1) = 0; chance2 = fittedpars; chance2(17:31) = 0;
%         [~, probs_ll] = getLogLik_step(truepars, data(s), false);
%         [~, probs_ll_fitted] = getLogLik_step(fittedpars, data(s), false);
%         [~, probs_choice] = makeChoice(truepars, data(s).options, data(s).avail_atts);
%         for i = 1:length(probs_ll), prob_real(i) = probs_ll(i,data(s).choices(i)); end
%         for i = 1:length(probs_ll), prob_real_fitted(i) = probs_ll_fitted(i,data(s).choices(i)); end

        [x,logpost] = ga(WAD_post, numParams,[],[],[],[],vertcat(param_struct_step.lb),...
            vertcat(param_struct_step.ub),[],(numParams-numAtts+1):numParams);

%         checkname = ['check_WAD_' num2str(s) '.mat'];
%         surrogateopts = optimoptions('surrogateopt','CheckpointFile',checkname,'MaxFunctionEvaluations',200);
%         disp(['Running initial for WAD_step for subject ', num2str(s)]);
%         logpost = 0;
%         [~,logpost2] = surrogateopt(WAD_post, lbs1, ubs1, (numParams-numAtts+1):numParams,surrogateopts);
%         num_reruns = 1;
%         while abs(logpost2 - logpost) > .01 && surrogateopts.MaxFunctionEvaluations < 1000
%             surrogateopts.MaxFunctionEvaluations = surrogateopts.MaxFunctionEvaluations + 50;
%             logpost = logpost2;
%             disp(['Running rerun #', num2str(num_reruns), ' for WAD_step for subject ', num2str(s), ', current ll = ', num2str(logpost)]);
%             [x,logpost2] = surrogateopt(checkname,surrogateopts);
%             num_reruns = num_reruns + 1;
%         end
%         logpost = logpost2;

        best_fit_params_WAD(s,:) = x;
        loglik_WAD(s) = getLogLik_step(x, data(s), false);
        logpost_WAD(s) = -logpost;
        BICs_WAD(s) = numParams*log(data(s).N) - 2*loglik_WAD(s);
        AICs_WAD(s) = numParams*2 - 2*loglik_WAD(s);

        disp(['Completed optimization for WAD_step for subject ', num2str(s)]);
    end

    % WP
    if any(strcmp(which_models, 'WP'))
        disp(['Fitting WP_step for subject ', num2str(s)]);
        WP_post = @(x) -post_WAD(x,data(s),true,param_struct_step,numAtts);
        [x,logpost] = ga(WP_post, numParams,[],[],[],[],lbs1,ubs1,[],(numParams-numAtts+1):numParams);
        best_fit_params_WP(s,:) = x;
        loglik_WP(s) = getLogLik(x, data(s), true);
        logpost_WP(s) = -logpost;
        BICs_WP(s) = numParams*log(data(s).N) - 2*loglik_WP(s);
        AICs_WP(s) = numParams*2 - 2*loglik_WP(s);
    end

    % EW
    if any(strcmp(which_models, 'EW'))
        disp(['Fitting EW_step for subject ', num2str(s)]);
        EW_post = @(x) -(getLogLik(x,data(s),false) + param_struct_step(1).logpdf(x(1)) + log((1/3) ^ numAtts) + log((1/(numAtts+1) ^ numAtts)));
        [x,logpost] = ga(EW_post, numParams,[],[],[],[],lbs2,ubs2,[],2:numParams);
        best_fit_params_EW(s,:) = x;
        loglik_EW(s) = getLogLik(x, data(s),false);
        logpost_EW(s) = -logpost;
        BICs_EW(s) = numParams*log(data(s).N) - 2*loglik_EW(s);
        AICs_EW(s) = numParams*2 - 2*loglik_EW(s);
    end

    % TAL
    if any(strcmp(which_models, 'TAL'))
        disp(['Fitting TAL_step for subject ', num2str(s)]);
        TAL_post = @(x) -(getLogLik(x,data(s),true) + param_struct_step(1).logpdf(x(1)) + log((1/3) ^ numAtts) + log((1/(numAtts+1) ^ numAtts)));

        [x,logpost] = ga(TAL_post, numParams,[],[],[],[],lbs2,ubs2,[],2:numParams);
        best_fit_params_TAL(s,:) = x;
        loglik_TAL(s) = getLogLik(x, data(s),true);
        logpost_TAL(s) = -logpost;
        BICs_TAL(s) = numParams*log(data(s).N) - 2*loglik_TAL(s);
        AICs_TAL(s) = numParams*2 - 2*loglik_TAL(s);
    end
end

results_template.K = numParams;
results_template.S = numSubj;
results_template.param = param_struct_step;
results_template.H = hessians;

results_WADs = results_template;
results_WADs.logpost = logpost_WAD';
results_WADs.loglik = loglik_WAD;
results_WADs.x = best_fit_params_WAD;
results_WADs.bic = BICs_WAD;
results_WADs.aic = AICs_WAD;
results_WADs.likfun = @(x,d) getLogLik_step(x, d, false);

results_WPs = results_template;
results_WPs.logpost = logpost_WP';
results_WPs.loglik = loglik_WP;
results_WPs.x = best_fit_params_WP;
results_WPs.bic = BICs_WP;
results_WPs.aic = AICs_WP;

results_EWs = results_template;
results_EWs.logpost = logpost_EW';
results_EWs.loglik = loglik_EW;
results_EWs.x = best_fit_params_EW;
results_EWs.bic = BICs_EW;
results_EWs.aic = AICs_EW;

results_TALs = results_template;
results_TALs.logpost = logpost_TAL';
results_TALs.loglik = loglik_TAL;
results_TALs.x = best_fit_params_TAL;
results_TALs.bic = BICs_TAL;
results_TALs.aic = AICs_TAL;

end