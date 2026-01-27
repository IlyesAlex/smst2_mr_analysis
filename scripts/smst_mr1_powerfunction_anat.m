%% smst_mr1_powerfunction_export_NOFUNCS.m
% Baseline-correct (AllFirst=0) per Subj×ROI, fit on ROI group means.
% Exports:
%   matlab_roi_means.csv
%   matlab_fit_table_long.csv
%   matlab_curve_predictions.csv

clear; clc;

%% 0) Load & normalize
inFile = "transfer_raw_matlab_anat.csv";
raw = readtable(inFile);

raw.Subj = string(raw.Subj);
raw.ROI  = string(raw.ROI);
raw.Cond = strtrim(string(raw.Cond));
raw.Beta = double(raw.Beta);

baselineCond = "AllFirst";
fitConds = ["ExactRepeat","CloseRepeat","DistantRepeat"];

%% 1) Baseline per Subj×ROI from AllFirst
baseRows = raw(raw.Cond == baselineCond, :);
if isempty(baseRows)
    error('No baseline rows with Cond == "%s". Check condition names.', baselineCond);
end

[gBase, baseSubj, baseROI] = findgroups(baseRows.Subj, baseRows.ROI);
baselineBeta = splitapply(@mean, baseRows.Beta, gBase);

baselineTbl = table(baseSubj, baseROI, baselineBeta, ...
    'VariableNames', {'Subj','ROI','baselineBeta'});

allWithBase = outerjoin(raw, baselineTbl, 'Keys', {'Subj','ROI'}, 'MergeKeys', true);
allWithBase.deltaBeta = allWithBase.Beta - allWithBase.baselineBeta;

%% 2) Keep fit conditions + map to x=1..3
fitRows = allWithBase(ismember(allWithBase.Cond, fitConds), :);

fitRows.x = NaN(height(fitRows),1);
fitRows.x(strcmpi(fitRows.Cond, "ExactRepeat"))   = 1;
fitRows.x(strcmpi(fitRows.Cond, "CloseRepeat"))   = 2;
fitRows.x(strcmpi(fitRows.Cond, "DistantRepeat")) = 3;

fitRows = fitRows(isfinite(fitRows.deltaBeta) & ~isnan(fitRows.x), :);

%% 3) ROI means ± SEM across subjects
[gMean, roiKey, xKey] = findgroups(fitRows.ROI, fitRows.x);

meanDelta = splitapply(@mean, fitRows.deltaBeta, gMean);
sdDelta   = splitapply(@std,  fitRows.deltaBeta, gMean);
nObs      = splitapply(@numel,fitRows.deltaBeta, gMean);
semDelta  = sdDelta ./ sqrt(nObs);

roiMeans = table(roiKey, xKey, meanDelta, semDelta, nObs, ...
    'VariableNames', {'ROI','x','mean_delta','sem_delta','n'});

roiMeans.Cond = strings(height(roiMeans),1);
roiMeans.Cond(roiMeans.x==1) = "ExactRepeat";
roiMeans.Cond(roiMeans.x==2) = "CloseRepeat";
roiMeans.Cond(roiMeans.x==3) = "DistantRepeat";

roiList = unique(roiMeans.ROI);

%% 4) Fit per ROI + export tables
ftPow = fittype('a*x^b', 'independent','x', 'coefficients',{'a','b'});
powOpts = fitoptions('Method','NonlinearLeastSquares');
powOpts.Display = 'Off';
try
    powOpts.Algorithm = 'Trust-Region';
catch
end

% Pre-allocate empty output tables
fitLong = table('Size',[0 9], ...
    'VariableTypes', {'string','string','string','double','double','double','double','double','string'}, ...
    'VariableNames', {'ROI','FitType','Model','AdjR2','R2','RMSE','SSE','DFE','Winner'});

curvePred = table('Size',[0 5], ...
    'VariableTypes', {'string','string','double','double','logical'}, ...
    'VariableNames', {'ROI','model','x','y','is_winner'});

xGrid = linspace(1,3,200)';

for i = 1:numel(roiList)
    roi = roiList(i);
    dd  = sortrows(roiMeans(roiMeans.ROI==roi, :), 'x');

    if height(dd) ~= 3
        warning("ROI %s has %d points (expected 3). Skipping.", roi, height(dd));
        continue;
    end

    x = dd.x;
    y = dd.mean_delta;

    % ----------------------------
    % Linear: polyfit
    % ----------------------------
    pLin = polyfit(x, y, 1);            % y = a*x + b
    aLin = pLin(1);
    bLin = pLin(2);
    yHatLin = polyval(pLin, x);

    % Fit stats (k=2)
    n = numel(y);
    sse = sum((y - yHatLin).^2);
    sst = sum((y - mean(y)).^2);
    if sst == 0
        r2 = NaN;
    else
        r2 = 1 - sse/sst;
    end
    dfe = n - 2;
    if isnan(r2) || dfe <= 0
        adjr2 = NaN; rmse = NaN;
    else
        adjr2 = 1 - (sse/dfe) / (sst/(n-1));
        rmse  = sqrt(sse/dfe);
    end

    linAdjR2 = adjr2; linR2 = r2; linRMSE = rmse; linSSE = sse; linDFE = dfe;
    linStr = sprintf("f(x)=%.4g*x + %.4g", aLin, bLin);

    % ----------------------------
    % Power2: fit a*x^b
    % ----------------------------
    aPow = NaN; bPow = NaN;
    powAdjR2 = NaN; powR2 = NaN; powRMSE = NaN; powSSE = NaN; powDFE = n - 2;
    powStr = "NA";

    try
        powFit = fit(x, y, ftPow, powOpts);
        aPow = powFit.a; bPow = powFit.b;
        yHatPow = aPow .* (x .^ bPow);

        sse = sum((y - yHatPow).^2);
        sst = sum((y - mean(y)).^2);
        if sst == 0
            r2 = NaN;
        else
            r2 = 1 - sse/sst;
        end
        dfe = n - 2;
        if isnan(r2) || dfe <= 0
            adjr2 = NaN; rmse = NaN;
        else
            adjr2 = 1 - (sse/dfe) / (sst/(n-1));
            rmse  = sqrt(sse/dfe);
        end

        powAdjR2 = adjr2; powR2 = r2; powRMSE = rmse; powSSE = sse; powDFE = dfe;
        powStr = sprintf("f(x)=%.4g*x^(%.4g)", aPow, bPow);
    catch ME
        warning("Power fit failed for ROI=%s (%s)", roi, ME.message);
    end

    % Winner by adjusted R^2
    if ~isnan(powAdjR2) && powAdjR2 > linAdjR2
        winner = "Power2";
    else
        winner = "Linear";
    end

    % Append fitLong (2 rows per ROI)
    fitLong = [fitLong; ...
        {roi, "Linear", linStr, linAdjR2, linR2, linRMSE, linSSE, linDFE, winner}; ...
        {roi, "Power2", powStr, powAdjR2, powR2, powRMSE, powSSE, powDFE, winner} ...
    ];

    % Dense predictions for plotting in R
    yGridLin = aLin .* xGrid + bLin;
    yGridPow = aPow .* (xGrid .^ bPow);

    curvePred = [curvePred; ...
        table(repmat(roi, numel(xGrid), 1), repmat("Linear", numel(xGrid), 1), xGrid, yGridLin, ...
              repmat(winner=="Linear", numel(xGrid), 1), ...
              'VariableNames', {'ROI','model','x','y','is_winner'}); ...
        table(repmat(roi, numel(xGrid), 1), repmat("Power2", numel(xGrid), 1), xGrid, yGridPow, ...
              repmat(winner=="Power2", numel(xGrid), 1), ...
              'VariableNames', {'ROI','model','x','y','is_winner'}) ...
    ];
end

disp(fitLong);

%% 5) Export
writetable(roiMeans, "matlab_roi_means_anat.csv");
writetable(fitLong,  "matlab_fit_table_long_anat.csv");
writetable(curvePred,"matlab_curve_predictions_anat.csv");

fprintf("\nExported:\n  matlab_roi_means_anat.csv\n  matlab_fit_table_long_anat.csv\n  matlab_curve_predictions_anat.csv\n");
