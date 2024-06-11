%% Script to simulate SSRN model response to pulsed input,
% and its modification by thalamic and cortical input
% adding basis functions of input from BG to thalamus
%
% Refactored: 3/3/2022
% 
% Mark Humphries 

clearvars; close all

% saving data
if ispc
    addpath('C:\Users\lpzmdh\Documents\GitHub\WriteSimulationData\');
else
    addpath('/Users/mqbssmhg/Documents/GitHub/WriteSimulationData/');
end

%% parameters for the RNN - overwrite SOC specific stuff with random network
run RNN_SOC_BG_TC_common_parameters

% load network
% networkID = 'Random_Network_20220210T104722';  % random 50:50 E:I
% networkID = 'Random_Network_20240610T160119_80_20'; % random 80:20 E:I
networkID = 'Random_Network_20240611T111018_80_20'; % another random 80:20 E:I
% networkID = 'SOC_Network_20240610T155910_80_20'; % SOC optimised 80:20 E:I

load(['Network building/' networkID]);

% set up output function
pars.output_type = 'power';
pars.output_arg1 = 0.1;            % scaling of output
pars.output_arg2 = 2;               % exponent of power function


if ispc
    savepath = 'C:\Users\lpzmdh\Dropbox\Simulations\Results\BG Mean-Field Model\Individual Sims\Basis functions\';
else
    savepath = '/Users/mqbssmhg/Dropbox/Simulations/Results/BG Mean-Field Model/Individual Sims/Basis functions/';
end


%% parameters for the TC and BG input to the RNN

pars.Pinput = 0.1;  % proportion of RNN units getting thalamic input
pars.ThalamicBaseline = pars.range(2);  % baseline output of thalamic neurons, including its weight

pars.BG_inputs = 5;  % number of BG inputs to thalamus
pars.BG_strength_steps = 5;     % N_s: number of discrete steps in BG weight distribution
pars.BG_baseline = pars.range(2) / 5;  % baseline of BG output, including its weight
pars.BG_deviation = pars.BG_baseline; % max deviation from baseline (assume symmetric for now)
pars.n_BG_tests = 10; % number of BG combinations to build and test

%% initialise TC inputs to cortex

% sample initial input vector to all RNN units
I_initial = pars.range(1) + (pars.range(2) - pars.range(1)) * rand(Net.N,1); 

% choose RNN units getting input
pars.ixThalamicInput = randi(Net.N,[round(Net.N*pars.Pinput),1]); % any unit: excite or inhibit

%% define BG output assignments to thalamus
number_of_Thalamus_neurons = numel(pars.ixThalamicInput);
% create template basis function, of width 2N_s - 1
template = zeros(1,number_of_Thalamus_neurons);
template(1:2*pars.BG_strength_steps - 1) = 1 - abs((1:2*pars.BG_strength_steps-1) - pars.BG_strength_steps) / pars.BG_strength_steps; % symmetric discrete basis function, peak of 1 

% create matrix of basis functions, one row per BG output neuron, on
% 1D-ring
for iB = 1:pars.BG_inputs
    pars.BG_basis_fcn_Weights(iB,:) = circshift(template,(iB-1)*(pars.BG_strength_steps-1));
end

% define BG output changes from baseline for each test
pars.BG_range = [pars.BG_baseline - pars.BG_deviation pars.BG_baseline + pars.BG_deviation];
pars.BG_inputs_to_thalamus = pars.BG_range(1) + (pars.BG_range(2) - pars.BG_range(1)) * rand(pars.BG_inputs,pars.n_BG_tests); % creates one vector of changes per column

%% create thalamic outputs, and update cortical input vector
pars.InputVectors = zeros(Net.N,pars.n_BG_tests); % storage: one column per input vector
    
for iInputs = 1:pars.n_BG_tests
    % compute total basal ganglia input to each thalamic neuron: multiply
    % row vector of BG output on this test with matrix of basis functions;
    % model: a'D = f
    BG_input_per_thalamic_neuron = pars.BG_inputs_to_thalamus(:,iInputs)' * pars.BG_basis_fcn_Weights; 
    % subtract from baseline the total BG input to each thalamic neuron
    ThalamicOutput = pars.ThalamicBaseline - BG_input_per_thalamic_neuron;
    % keyboard
    % assign new thalamic output to input vector
    pars.InputVectors(:,iInputs) = I_initial;
    pars.InputVectors(pars.ixThalamicInput,iInputs) = ThalamicOutput;
end

%% run network
tSteps = round(pars.Tmax ./ pars.dt);   % how many time-steps

% time duration of input too
tInput = round(pars.duration ./pars.dt);

% storage
Data.rStore = zeros(pars.n_BG_tests*tSteps,Net.N);  % unit outputs
Data.temporalPCs = zeros(pars.n_BG_tests,tSteps,pars.nPCs); % PC projections
    
for iBG = 1:pars.n_BG_tests
    
    % run simulation with that vector
    a = zeros(Net.N,tSteps);
    r = zeros(Net.N,tSteps);

    for iTime = 2:tSteps

        % background input
        I = zeros(Net.N,1) + pars.I_background;

        % add stepped input
        if iTime > tInput(1) && iTime <= tInput(2)
            I = I + pars.InputVectors(:,iBG);
        end
        % update activity
        a(:,iTime) = a(:,iTime-1) + pars.dt * (-a(:,iTime-1) + I + optimW*r(:,iTime-1)) ./ pars.ts;

        % update output
        r(:,iTime) = neuron_output(a(:,iTime),pars.output_type,pars.output_arg1,pars.output_arg2);

    end

    % store outputs, concatenated
    Data.rStore(1+(iBG-1)*tSteps:iBG*tSteps,:) = r(:,1:tSteps)';  % columns are neurons
        
end

% compute PCA on all concatenated outputs
[Vectors,temporalPCs,EigVals] = pca(Data.rStore);   % do PCA on all neuron outputs


% store top N PCs, for each BG input, stored as a matrix
for iBG = 1:pars.n_BG_tests
    Data.temporalPCs(iBG,:,:) = temporalPCs(1+(iBG-1)*tSteps:iBG*tSteps,1:pars.nPCs);
end

%% compute Euclidean distance between pairs of trajectories 
% compare to distance between input vectors and BG vectors
Data.maxdistance = zeros(pars.n_BG_tests);
for iBG = 1:pars.n_BG_tests
    for jBG = iBG+1:pars.n_BG_tests
%             tic
%             for iTime = 1:tSteps
%                 % distance = sqrt(sum((a-b)^2))
%                 diffs = squeeze(bsxfun(@minus,Data.temporalPCs(iBatch,iAngle,iTime,:),Data.temporalPCs(iBatch,jAngle,iTime,:)));
%                 distance(iAngle,jAngle,iTime) = sqrt(sum(diffs.^2));
%             end
%            toc

        % compute Euclidean distance between trajectories
        for iPC = 1:pars.nPCs
            % for each PC compute the difference between the
            % time-points at each angle
            diffsPCs(iPC,:) = squeeze(bsxfun(@minus,Data.temporalPCs(iBG,:,iPC),Data.temporalPCs(jBG,:,iPC)));
        end
        % compute Euclidean distance from those differences
        distancePCs(iBG,jBG,:) = sqrt(sum(diffsPCs.^2));

        % save max distance between trajectories for analysis
        Data.maxdistance(iBG,jBG) = max(distancePCs(iBG,jBG,:));
        Data.maxdistance(jBG,iBG) = Data.maxdistance(iBG,jBG);
        
    end
end

% compute Euclidean distance between all input vectors - pdist
% computes between rows
Data.InputVectorDistances = squareform(pdist(pars.InputVectors'));
Data.BGOutputDistances = squareform(pdist(pars.BG_inputs_to_thalamus'));

%% save data
F = makeRecordOfSimulationData('RandomNetwork_BG_BasisFcn_',savepath,Data,pars,tSteps,tInput,networkID);

%% plot checks
rate_color_map = brewermap(15,'Greys');
figure
imagesc(Data.rStore')
colormap(rate_color_map)
colorbar
xlabel('Time')
ylabel('Neuron')


PCcmap = brewermap(pars.n_BG_tests,'PuOr');
% plot PCs on same figure
nPlots = 3;
figure
for iP = 1:nPlots
    plot3(squeeze(Data.temporalPCs(iP,:,1)),squeeze(Data.temporalPCs(iP,:,2)),squeeze(Data.temporalPCs(iP,:,3)),'Color',PCcmap(iP,:));
    hold on
end

% plot distances between inputs and outputs
figure
for iBG = 1:pars.n_BG_tests
    plot(Data.InputVectorDistances(iBG,:),Data.maxdistance(iBG,:),'o'); hold on
end
xlabel('Distance between input vectors')
ylabel('Max distance between projections')

figure
for iBG = 1:pars.n_BG_tests
    plot(Data.BGOutputDistances(iBG,:),Data.maxdistance(iBG,:),'o'); hold on
end
xlabel('Distance between BG output vectors')
ylabel('Max distance between projections')








