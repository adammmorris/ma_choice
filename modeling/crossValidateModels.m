function [cv_WAD, cv_WP, cv_EW, cv_TAL, cv_LEX] = crossValidateModels(data, results_all)
cv_WAD = mfit_predict(data, results_all(1));
cv_WP = mfit_predict(data, results_all(2));
cv_EW = mfit_predict(data, results_all(3));
cv_TAL = mfit_predict(data, results_all(4));
cv_LEX = mfit_predict(data, results_all(5));

end