function [SNr,mdl_upper,mdl_lower,mdl_weighted_upper,mdl_weighted_lower] = plot_and_fit(Connections,InjectionVolume,volume_of_SNr,Proportion_in_SNr)

% get data for estimating scaling
DVupper = [Connections(:).totalNeurons];     
DVlower = [Connections(:).totalLowerBound];

IV = InjectionVolume; % use raw injection volume if looking at experiments outside SNr
xText = 'Injection volume (mm^3)';


% do regression on volume versus total number of neurons - Stats toolbox needed for linear
% model functions (fitlm; predict)
mdl_upper = fitlm(IV',DVupper);
mdl_lower = fitlm(IV',DVlower);

weights = Proportion_in_SNr ./ max(Proportion_in_SNr);        % weighted by proportion of injection in SNr
mdl_weighted_upper = fitlm(IV',DVupper,'Weights',weights);
mdl_weighted_lower = fitlm(IV',DVlower,'Weights',weights);

% create regression lines to plot
evaluate_between = [min(InjectionVolume)-0.05:0.02:volume_of_SNr+0.05]';  % start and end points to plot regression line
[UpperBound_fit,UpperBoundCI] = predict(mdl_upper,evaluate_between);
[UpperBound_weighted_fit,UpperBound_weighted_CI] = predict(mdl_weighted_upper,evaluate_between);

[LowerBound_fit,LowerBoundCI] = predict(mdl_lower,evaluate_between);
[LowerBound_weighted_fit,LowerBound_weighted_CI] = predict(mdl_weighted_lower,evaluate_between);

% look up SNr-scale estimate
[SNr.upper_estimate,SNr.upper_estimate_CI] = predict(mdl_upper,volume_of_SNr);
[SNr.weighted_upper_estimate,SNr.weighted_upper_estimate_CI] = predict(mdl_weighted_upper,volume_of_SNr);

[SNr.lower_estimate,SNr.lower_estimate_CI] = predict(mdl_lower,volume_of_SNr);
[SNr.weighted_lower_estimate,SNr.weighted_lower_estimate_CI] = predict(mdl_weighted_lower,volume_of_SNr);

% plot numbers
figure
subplot(221),plot(IV,DVlower,'k.'); hold on
plot(evaluate_between,LowerBound_fit,'Color',[0.8 0.3 0.3])
plot(evaluate_between,LowerBoundCI,'Color',[0.8 0.7 0.7])
axisH = gca;
line([volume_of_SNr volume_of_SNr],[axisH.YLim(1) SNr.lower_estimate])
line([axisH.XLim(1) volume_of_SNr],[SNr.lower_estimate SNr.lower_estimate])

ylabel('Targets [density-corrected]')
xlabel(xText)
title('Lower Bounds')

subplot(223),scatter(IV,DVlower,weights*50); hold on
plot(evaluate_between,LowerBound_weighted_fit,'Color',[0.8 0.3 0.3])
plot(evaluate_between,LowerBound_weighted_CI,'Color',[0.8 0.7 0.7])
axisH = gca;
line([volume_of_SNr volume_of_SNr],[axisH.YLim(1) SNr.weighted_lower_estimate])
line([axisH.XLim(1) volume_of_SNr],[SNr.weighted_lower_estimate SNr.weighted_lower_estimate])
%ylabel('Ratio target:SNr [density-corrected]')
xlabel(xText)

% plot upper bounds
subplot(222)
plot(IV,DVupper,'k.'); hold on
plot(evaluate_between,UpperBound_fit,'Color',[0.8 0.3 0.3])
plot(evaluate_between,UpperBoundCI,'Color',[0.8 0.7 0.7])
axisH = gca;
line([volume_of_SNr volume_of_SNr],[axisH.YLim(1) SNr.upper_estimate])
line([axisH.XLim(1) volume_of_SNr],[SNr.upper_estimate SNr.upper_estimate])

ylabel('Total target neurons')
xlabel(xText)
title('Upper Bound')

subplot(224),scatter(IV,DVupper,weights*50); hold on
plot(evaluate_between,UpperBound_weighted_fit,'Color',[0.8 0.3 0.3])
plot(evaluate_between,UpperBound_weighted_CI,'Color',[0.8 0.7 0.7])
axisH = gca;
line([volume_of_SNr volume_of_SNr],[axisH.YLim(1) SNr.weighted_upper_estimate])
line([axisH.XLim(1) volume_of_SNr],[SNr.weighted_upper_estimate SNr.weighted_upper_estimate])
xlabel(xText)