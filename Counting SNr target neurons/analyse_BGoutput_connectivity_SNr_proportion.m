% script to assess BG output target counts in all Allen brain experiments
% By looking at how target neuron count scales with approximate number of
% SNr neurons injected
% Can only do this for injections largely confined to SNr (otherwise cannot
% correct)
% 
% Mark Humphries 9/5/2024

clearvars; close all;

% which experiments to use
expts = [100141993, 175263063, 299895444]; % all with >90% injection in SNr

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

%% specify density threshold to assess
density_threshold = 0:0.005:0.05;

%% assess them...

resultsTableLower = zeros(numel(density_threshold),numel(expts));
resultsTableUpper = zeros(numel(density_threshold),numel(expts));

for ixParameters = 1:numel(density_threshold)
    Parameters(ixParameters).density_threshold = density_threshold(ixParameters);
    % compute connectivity
    for iExpt = 1:numel(expts)
        expt_filename = ['Experiment_' num2str(expts(iExpt)) '.xml'];
        Connections(iExpt).Data = count_target_neurons(expt_filename,Parameters(ixParameters).density_threshold);
        Connections(iExpt).totalNeurons = sum(Connections(iExpt).Data.NeuronsInTarget,"omitnan");
        Connections(iExpt).totalLowerBound = sum(Connections(iExpt).Data.LowerBound,"omitnan");

        % compute ratios for this experiment
        estimated_neurons_hit_in_SNr = neurons_in_SNr * InjectionVolume(iExpt) / volume_of_SNr; % estimate of how many SNr neurons it hit
        Connections(iExpt).ExpansionUpperBound = Connections(iExpt).totalNeurons ./ estimated_neurons_hit_in_SNr;
        Connections(iExpt).ExpansionLowerBound = Connections(iExpt).totalLowerBound  ./ estimated_neurons_hit_in_SNr;

        % count number of missing target structures from BBP database
        Connections(iExpt).totalTargets = height(Connections(iExpt).Data);
        Connections(iExpt).missingTargets = sum(isnan(Connections(iExpt).Data.NeuronsInTarget));

    end
    Parameters(ixParameters).Connections = Connections;
    resultsTableUpper(ixParameters,:) = [Connections(:).ExpansionUpperBound];
    resultsTableLower(ixParameters,:) = [Connections(:).ExpansionLowerBound];
end

TableResultsUpperBound = array2table(resultsTableUpper,"VariableNames",string(expts));
TableResultsUpperBound.density_threshold = density_threshold';
TableResultsLowerBound = array2table(resultsTableLower,"VariableNames",string(expts));
TableResultsLowerBound.density_threshold = density_threshold';

%% plot density threshold versus ratio

figure
subplot(211)
for ixParameters = 1:numel(Parameters)
    plot(Parameters(ixParameters).density_threshold,[Parameters(ixParameters).Connections(:).ExpansionLowerBound],'.'); hold on
end
set(gca,'XLim',[min([Parameters(:).density_threshold])-0.005 ,max([Parameters(:).density_threshold])+0.01])
xlabel('Density threshold')
ylabel('Lower target:SNr ratio')
subplot(212)
for ixParameters = 1:numel(Parameters)
    plot(Parameters(ixParameters).density_threshold,[Parameters(ixParameters).Connections(:).ExpansionUpperBound],'.'); hold on
end
xlabel('Density threshold')
ylabel('Upper target:SNr ratio')

%% find asymptotic values
eta = 1;
for iExpt = 1:numel(expts)
    [upper_bound(iExpt).data_asymptote,upper_bound(iExpt).fit_asymptote,upper_bound(iExpt).fitted_model] = find_asymptotic_ratio(TableResultsUpperBound.density_threshold,TableResultsUpperBound.(string(expts(iExpt))),eta);
    [lower_bound(iExpt).data_asymptote,lower_bound(iExpt).fit_asymptote,lower_bound(iExpt).fitted_model] = find_asymptotic_ratio(TableResultsLowerBound.density_threshold,TableResultsLowerBound.(string(expts(iExpt))),eta);
end


% final answer
lower_bound_data = mean([lower_bound(:).data_asymptote])
upper_bound_data = mean([upper_bound(:).data_asymptote])

save BG_output_connections_SNr_proportion expts Parameters TableResultsLowerBound TableResultsUpperBound upper_bound lower_bound
