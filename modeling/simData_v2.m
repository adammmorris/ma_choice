%% setup
clear
addpath("mfit/");

numAgents = 100;
numAtts = 15;
numParams = numAtts + 1;
numChoices = 75;
numValuesPerAtt = 5;

%% generate options
options = randi([1 numValuesPerAtt], numAtts, 2, numChoices, numAgents) / numValuesPerAtt;
options_signed = signOptions(options);

%% generate agent parameters
weight_params = [0 1];
weights = normrnd(weight_params(1),weight_params(2),numAgents,numAtts);
[~, max_weight_ind] = max(abs(weights),[],2);
weights_signed = sign(weights);

% I know there's a vectorized way to do this, but can't figure it out rn
weights_lex = zeros(numAgents, numAtts);
for i = 1:numAgents
    weights_lex(i, max_weight_ind(i)) = sign(weights(i, max_weight_ind(i)));
end

gamma_bounds = [1 5];
inv_temp = gamrnd(gamma_bounds(1),gamma_bounds(2),numAgents,1);

struct_template.N = numChoices;
struct_template.options = options;

for agent = 1:numAgents
    data_WAD(agent) = struct_template;
    data_WP(agent) = struct_template;
    data_EW(agent) = struct_template;
    data_TAL(agent) = struct_template;
    data_LEX(agent) = struct_template;
end

params_WAD = zeros(numAgents, numParams);
params_WP = zeros(numAgents, numParams);
params_EW = zeros(numAgents, numParams);
params_TAL = zeros(numAgents, numParams);
params_LEX = zeros(numAgents, numParams);

%% simulate data
for agent = 1:numAgents

    data_WAD(agent).options = options(:,:,:,agent);
    data_WAD(agent).params = [inv_temp(agent) weights(agent,:)];
    params_WAD(agent, :) = data_WAD(agent).params;
    data_WAD(agent).choices = makeChoice_WAD(data_WAD(agent).params, data_WAD(agent).options);

    data_WP(agent).options = options(:,:,:,agent);
    data_WP(agent).params = [inv_temp(agent) weights(agent,:)];
    params_WP(agent, :) = data_WP(agent).params;
    data_WP(agent).choices = makeChoice_WAD(data_WP(agent).params, signOptions(data_WP(agent).options));

    data_EW(agent).options = options(:,:,:,agent);
    data_EW(agent).params = [inv_temp(agent) weights_signed(agent,:)];
    params_EW(agent, :) = data_EW(agent).params;
    data_EW(agent).choices = makeChoice_WAD(data_EW(agent).params, data_EW(agent).options);

    data_TAL(agent).options = options(:,:,:,agent);
    data_TAL(agent).params = [inv_temp(agent) weights_signed(agent,:)];
    params_TAL(agent, :) = data_TAL(agent).params;
    data_TAL(agent).choices = makeChoice_WAD(data_TAL(agent).params, signOptions(data_TAL(agent).options));

    data_LEX(agent).options = options(:,:,:,agent);
    data_LEX(agent).params = [inv_temp(agent) weights_lex(agent,:)];
    params_LEX(agent, :) = data_LEX(agent).params;
    data_LEX(agent).choices = makeChoice_WAD(data_LEX(agent).params, data_LEX(agent).options);
end

save('simdata.mat');

%% fit models
param(1).name = 'inverse temperature';
param(1).logpdf = @(x) sum(log(gampdf(x,gamma_bounds(1),gamma_bounds(2))));  % log density function for prior
param(1).lb = 0;    % lower bound
param(1).ub = 50;   % upper bound

for i = 1:numAtts
    param(i+1).name = strcat('weight',string(i));
    param(i+1).logpdf = @(x) sum(log(normpdf(x,weight_params(1),weight_params(2))));  % log density function for prior
    param(i+1).lb = -5;    % lower bound
    param(i+1).ub = 5;   % upper bound
end

nstarts = 5;

[results_WAD_WAD, results_WAD_WP, results_WAD_EW, results_WAD_TAL, results_WAD_LEX] = ...
    fitModels(param, data_WAD, nstarts, numAtts);
[results_WP_WAD, results_WP_WP, results_WP_EW, results_WP_TAL, results_WP_LEX] = ...
    fitModels(param, data_WP, nstarts, numAtts);
[results_EW_WAD, results_EW_WP, results_EW_EW, results_EW_TAL, results_EW_LEX] = ...
    fitModels(param, data_EW, nstarts, numAtts);
[results_TAL_WAD, results_TAL_WP, results_TAL_EW, results_TAL_TAL, results_TAL_LEX] = ...
    fitModels(param, data_TAL, nstarts, numAtts);
[results_LEX_WAD, results_LEX_WP, results_LEX_EW, results_LEX_TAL, results_LEX_LEX] = ...
    fitModels(param, data_LEX, nstarts, numAtts);

save('fitting.mat');

%% Test fits
cors_WAD = zeros(numParams,1);
for i = 1:numParams
    cors_WAD(i) = corr(params_WAD(:,i), results_WAD_WAD.x(:,i));
end
cors_WAD(1)
mean(cors_WAD(2:end))

cors_WP = zeros(numParams,1);
for i = 1:numParams
    cors_WP(i) = corr(params_WP(:,i), results_WP_WP.x(:,i));
end
cors_WP(1)
mean(cors_WP(2:end))

cors_EW = zeros(numParams,1);
for i = 1:numParams
    cors_EW(i) = corr(params_EW(:,i), results_EW_EW.x(:,i));
end
cors_EW(1)
mean(cors_EW(2:end))

cors_TAL = zeros(numParams,1);
for i = 1:numParams
    cors_TAL(i) = corr(params_TAL(:,i), results_TAL_TAL.x(:,i));
end
cors_TAL(1)
mean(cors_TAL(2:end))

cors_LEX = zeros(numParams,1);
for i = 1:numParams
    cors_LEX(i) = corr(params_LEX(:,i), results_LEX_LEX.x(:,i));
end
cors_LEX(1)
mean(cors_LEX(2:end))

bms_results_WAD = mfit_bms([results_WAD_WAD, results_WAD_WP, results_WAD_EW, results_WAD_TAL, results_WAD_LEX]);
bms_results_WAD.pxp
bms_results_WP = mfit_bms([results_WP_WAD, results_WP_WP, results_WP_EW, results_WP_TAL, results_WP_LEX]);
bms_results_WP.pxp
bms_results_EW = mfit_bms([results_EW_WAD, results_EW_WP, results_EW_EW, results_EW_TAL, results_EW_LEX]);
bms_results_EW.pxp
bms_results_TAL = mfit_bms([results_TAL_WAD, results_TAL_WP, results_TAL_EW, results_TAL_TAL, results_TAL_LEX]);
bms_results_TAL.pxp
bms_results_LEX = mfit_bms([results_LEX_WAD, results_LEX_WP, results_LEX_EW, results_LEX_TAL, results_LEX_LEX]);
bms_results_LEX.pxp