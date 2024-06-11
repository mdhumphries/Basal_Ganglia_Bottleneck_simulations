% script to assess BG output target counts in all Allen brain experiments
% Mark Humphries 23/4/2024

clearvars; close all;

% which experiments to use
expts = [100141993, 158914182, 175263063, 299895444, 478096249, 478097069]; % all with decent sized injection volumes

% get data for each experiment
ExperimentMetaData = readtable('AllenMouseProjectionAtlas_SNr_experiments');

for iExpt = 1:numel(expts)
    InjectionVolume(iExpt) = ExperimentMetaData.injection_volume(ExperimentMetaData.id == expts(iExpt));
    ProportionInSNr(iExpt) = ExperimentMetaData.FractionOfEachInjectionInSNr_myAddition_(ExperimentMetaData.id == expts(iExpt));
end

%% as reference, compute number of neurons in SNr and striatum
BBP_Cell_Atlas_Data = readtable('Rodarie2023_PLoSCB_data.xlsx');
BBP_Cell_Atlas_Data.NeuronCounts = BBP_Cell_Atlas_Data.Neuron_mm__3_ .* BBP_Cell_Atlas_Data.Volumes_mm_3_;

indexSNr = strcmp(BBP_Cell_Atlas_Data.BrainRegion,'Substantia nigra, reticular part'); 
neurons_in_SNr = BBP_Cell_Atlas_Data.NeuronCounts(indexSNr) / 2;
volume_of_SNr = BBP_Cell_Atlas_Data.Volumes_mm_3_(indexSNr) / 2;

indexDorsalStriatum = strcmp(BBP_Cell_Atlas_Data.BrainRegion,'Striatum dorsal region'); 
neurons_in_DorsalStriatum = BBP_Cell_Atlas_Data.NeuronCounts(indexDorsalStriatum) / 4;  
% why 4? 1/2 to get single hemisphere; then 1/2 to get D1-receptor population projecting to SNr

Compression = neurons_in_DorsalStriatum ./ neurons_in_SNr;

%% specify density to assess
density_threshold = 0:0.005:0.05;


%% assess them...
for ixParameters = 1:numel(density_threshold)
    Parameters(ixParameters).density_threshold = density_threshold(ixParameters);
    % compute connectivity
    for iExpt = 1:numel(expts)
        expt_filename = ['Experiment_' num2str(expts(iExpt)) '.xml'];
        Connections(iExpt).Data = count_target_neurons(expt_filename,Parameters(ixParameters).density_threshold);
        Connections(iExpt).totalNeurons = sum(Connections(iExpt).Data.NeuronsInTarget,"omitnan");
        Connections(iExpt).totalLowerBound = sum(Connections(iExpt).Data.LowerBound,"omitnan");

        % count number of missing target structures from BBP database
        Connections(iExpt).totalTargets = height(Connections(iExpt).Data);
        Connections(iExpt).missingTargets = sum(isnan(Connections(iExpt).Data.NeuronsInTarget));

    end
    Parameters(ixParameters).Connections = Connections;
    [Parameters(ixParameters).SNrEstimates,Parameters(ixParameters).mdl_upper,Parameters(ixParameters).mdl_lower,...
        Parameters(ixParameters).weighted_mdl_upper,Parameters(ixParameters).weighted_mdl_lower] = plot_and_fit(Connections,InjectionVolume,volume_of_SNr,ProportionInSNr);
   
    % compute ratios given scaling of target numbers with injection volume
    Parameters(ixParameters).ExpectedUpperBoundRatio = Parameters(ixParameters).SNrEstimates.upper_estimate / neurons_in_SNr;
    Parameters(ixParameters).CIUpperBoundRatio = Parameters(ixParameters).SNrEstimates.upper_estimate_CI ./ neurons_in_SNr;
    Parameters(ixParameters).ExpectedWeightedUpperBoundRatio = Parameters(ixParameters).SNrEstimates.weighted_upper_estimate ./ neurons_in_SNr;
    Parameters(ixParameters).CIWeightedUpperBoundRatio = Parameters(ixParameters).SNrEstimates.weighted_upper_estimate_CI ./ neurons_in_SNr;
    Parameters(ixParameters).ExpectedLowerBoundRatio = Parameters(ixParameters).SNrEstimates.lower_estimate ./ neurons_in_SNr;
    Parameters(ixParameters).CILowerBoundRatio = Parameters(ixParameters).SNrEstimates.lower_estimate_CI ./ neurons_in_SNr;
    Parameters(ixParameters).ExpectedWeightedLowerBoundRatio = Parameters(ixParameters).SNrEstimates.weighted_lower_estimate ./ neurons_in_SNr;
    Parameters(ixParameters).CIWeightedLowerBoundRatio = Parameters(ixParameters).SNrEstimates.weighted_lower_estimate_CI ./ neurons_in_SNr;
end


%% plot ratio vs density threshold

% specify raw and weighted colours
% add CIs in light colours of each

figure
subplot(121)
semilogy([Parameters(:).density_threshold],[Parameters(:).ExpectedUpperBoundRatio]); hold on
semilogy([Parameters(:).density_threshold],[Parameters(:).ExpectedWeightedUpperBoundRatio],'color',[0.7 0.7 0.7]);
% ylabel('SNr:target ratio')
xlabel('Density threshold')
title('Upper bound')

subplot(122)
plot([Parameters(:).density_threshold],[Parameters(:).ExpectedLowerBoundRatio]); hold on
plot([Parameters(:).density_threshold],[Parameters(:).ExpectedWeightedLowerBoundRatio],'color',[0.7 0.7 0.7]);
ylabel('SNr:target ratio')
xlabel('Density threshold')
title('lower bound')

%% find asymptote
eta = 1;
[upper_bound.data_asymptote,upper_bound.fit_asymptote,upper_bound.fitted_model] = find_asymptotic_ratio([Parameters(:).density_threshold],[Parameters(:).ExpectedUpperBoundRatio],eta);
[upper_bound.weighted_data_asymptote,upper_bound.weighted_fit_asymptote,upper_bound.weighted_fitted_model] = find_asymptotic_ratio([Parameters(:).density_threshold],[Parameters(:).ExpectedWeightedUpperBoundRatio],eta);
[lower_bound.data_asymptote,lower_bound.fit_asymptote,lower_bound.fitted_model] = find_asymptotic_ratio([Parameters(:).density_threshold],[Parameters(:).ExpectedLowerBoundRatio],eta);
[lower_bound.weighted_data_asymptote,lower_bound.weighted_fit_asymptote,lower_bound.weighted_fitted_model] = find_asymptotic_ratio([Parameters(:).density_threshold],[Parameters(:).ExpectedWeightedLowerBoundRatio],eta);

save BG_output_connections_scaling expts Parameters upper_bound lower_bound