%% setup
clear
datapath = '../data/v2/pilot2/';

addpath("mfit/");
load(strcat(datapath, 'imported_data.mat'));

param(1).name = 'inverse temperature';
param(1).logpdf = @(x) sum(log(gampdf(x,1,5)));  % log density function for prior
%param(1).logpdf = @(x) 0;
param(1).lb = 0;    % lower bound
param(1).ub = 50;   % upper bound

for i = 1:numAtts
    param(i+1).name = strcat('weight',string(i));
    param(i+1).logpdf = @(x) sum(log(normpdf(x,0,1)));  % log density function for prior
    %param(i+1).logpdf = @(x) 0;
    param(i+1).lb = -5;    % lower bound
    param(i+1).ub = 5;   % upper bound
end

nstarts = 5;

[results_WAD, results_WP, results_EW, results_TAL] = ...
    fitModels(param, data_real, nstarts, numAtts);

results_all = [results_WAD, results_WP, results_EW, results_TAL];
writematrix(results_WAD.x(:,2:end), strcat(datapath, 'fitted_empirical_weights_WAD.csv'))
writematrix(results_WP.x(:,2:end), strcat(datapath, 'fitted_empirical_weights_WP.csv'))
writematrix(results_EW.x(:,2:end), strcat(datapath, 'fitted_empirical_weights_EW.csv'))
writematrix(results_TAL.x(:,2:end), strcat(datapath, 'fitted_empirical_weights_TAL.csv'))

% do bayesian model selection
bms_results = mfit_bms(results_all, 1);
bms_results.pxp

[results_WAD.bic results_WP.bic results_EW.bic results_TAL.bic]

% do cross-validation

numFolds = 5;
numModels = 4;

for subj = 1:numSubj
    numChoices_cur = data_real(subj).N;
    numTrainTrials = ceil(numChoices_cur - numChoices_cur / numFolds);
    numTestTrials = floor(numChoices_cur / numFolds);
    choices_shuffled = 1:numChoices_cur;
    choices_shuffled = choices_shuffled(randperm(length(choices_shuffled)));

    for i = 1:numFolds
        test_trials = choices_shuffled((numTestTrials*(i-1)+1):(numTestTrials*i));
        train_trials = choices_shuffled;
        train_trials(test_trials) = [];

        traindata(subj).N = numTrainTrials;
        traindata(subj).options = data_real(subj).options(:,:,train_trials);
        traindata(subj).avail_atts = data_real(subj).avail_atts(train_trials,:);
        traindata(subj).choices = data_real(subj).choices(train_trials,:);

        testdata(subj).N = numTestTrials;
        testdata(subj).options = data_real(subj).options(:,:,test_trials);
        testdata(subj).avail_atts = data_real(subj).avail_atts(test_trials,:);
        testdata(subj).choices = data_real(subj).choices(test_trials,:);

        folds(i).traindata = traindata;
        folds(i).testdata = testdata;
    end
end

cv_results = zeros(numSubj, numModels, numFolds);
for i = 1:numFolds
    [results_WAD_train, results_WP_train, results_EW_train, results_TAL_train] = ...
        fitModels(param, folds(i).traindata, nstarts, numAtts);
    results_all_train = [results_WAD_train, results_WP_train, results_EW_train, results_TAL_train];
    [results_WAD_test, results_WP_test, results_EW_test, results_TAL_test] = ...
        crossValidateModels(folds(i).testdata, results_all_train);
    results_all_test = [results_WAD_test, results_WP_test, results_EW_test, results_TAL_test];

    cv_results(:,:,i) = results_all_test;
end

cv_results_avg = mean(cv_results,3);
mean(mean(cv_results,3))
median(mean(cv_results,3))

best_model = zeros(numSubj, 1);
for i = 1:numSubj
    [~, best_model(i)] = max(cv_results_avg(i,:));
end
graphCV(cv_results_avg, [1 3]);

hist(best_model)

save(strcat(datapath,'fit_empirical.mat'));

%% normalize
chance = log(.5 ^ 15);
cv_results_normalized = zeros(numSubj, numModels);
for i = 1:numSubj
    best = cv_results_avg(i,best_model(i));
    cv_results_normalized(i,:) = (cv_results_avg(i,:) - chance) / (best - chance);
end

writematrix(cv_results_normalized, strcat(datapath, 'cv_results.csv'))
