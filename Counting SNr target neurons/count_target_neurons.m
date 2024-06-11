function NeuronCountTable = count_target_neurons(experiment,density_threshold)

% CountTargetNeurons counting neurons in target structures of the SNr
% T = CountTargetNeurons(EXPT,THRESHOLD) analyses the Allen Mouse
% Brain Connectivity Atlas data in XML file EXPT. It selects target regions
% with a projection density *greater* than THRESHOLD (proportion
% between 0 and 1) in either hemisphere. 
% 
% It then selects the targets that appear in the non-overlapping list of
% brain regions from the Oh et al (2014) Nature paper
%
% These regions are then looked up in the Blue Brain Project cell count
% atlas (Rodrie et al 2023 PloS CB) to find their target numbers
% 
% A lower bound on connections is computed by multiplying these target
% bounds by the projection density
%
% 23/4/2024: initial version
% 01/5/2024: updated to use full target list from Oh et al
% 
% Mark Humphriess

AllenCCFData = readtable('Complete_Allen_Mouse_CCF_namespace.xlsx'); % load namespace for IDs of regions

%% load and process connectivity data
ConnectivityData = readtable(experiment);

potential_targets = find((ConnectivityData.("hemisphere-id") == 1 | ConnectivityData.("hemisphere-id") == 2)  & ...  % check targets in both hemispheres 
                        ConnectivityData.("is-injection") == "false" & ...            % check entry is not counting injection itself
                        ConnectivityData.("projection-density") > density_threshold);   % check that target passes threshold check
potential_TargetIDs = ConnectivityData.("structure-id")(potential_targets);

% filter by set of targets in Oh et al 2014 data
Allen_paper_regions = readtable('Oh2014_regions_for_connectome.xlsx');

[TargetIDs,index_into_potential_targets] = intersect(potential_TargetIDs,Allen_paper_regions.ID); % keep all TargetIDs in Allen_paper_regions list
targets_to_use = potential_targets(index_into_potential_targets);

% find names for all targets and their depth in the CCF tree
TargetName = cell(size(TargetIDs)); TargetAbbrev = TargetName; 
for ixID = 1:numel(TargetIDs) % map to names
    index_in_CCF_table = AllenCCFData.structureID == TargetIDs(ixID);
    TargetName{ixID} = AllenCCFData.fullStructureName{index_in_CCF_table};
    TargetAbbrev{ixID} = AllenCCFData.abbreviation{index_in_CCF_table};
end

%% count neurons in targets
BBP_Cell_Atlas_Data = readtable('Rodarie2023_PLoSCB_data.xlsx');
BBP_Cell_Atlas_Data.NeuronCounts = BBP_Cell_Atlas_Data.Neuron_mm__3_ .* BBP_Cell_Atlas_Data.Volumes_mm_3_;

NeuronsInTarget = zeros(size(TargetIDs));
for ixID = 1:numel(TargetIDs)
    index_in_BBP_Data = strcmp(BBP_Cell_Atlas_Data.BrainRegion,TargetName{ixID}); % find exact match to full name
    if sum(index_in_BBP_Data) == 1
        NeuronsInTarget(ixID) = BBP_Cell_Atlas_Data.NeuronCounts(index_in_BBP_Data) / 2;  % divide by 2 to get numbers in 1 hemisphere
    elseif sum(index_in_BBP_Data) == 0  % BBP data has no entry 
        NeuronsInTarget(ixID) = nan;
    else % something gone horribly wrong
        keyboard
    end
end

%% do lower bound: neurons * density
LowerBound = NeuronsInTarget .* ConnectivityData.("projection-density")(targets_to_use);

%% create table
NeuronCountTable = table(TargetIDs,TargetName,TargetAbbrev,NeuronsInTarget,LowerBound);
NeuronCountTable.Properties.VariableNames = {'TargetIDs','TargetName','TargetAbbrev','NeuronsInTarget','LowerBound'};




