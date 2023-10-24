%% script to create random balanced network with specified spectral radius 
% and save results
%
% Mark Humphries 10/2/2022
clearvars; close all

%% initialise
% parameterise network
Net.N = 200;
Net.p = 0.1;    % probability of connection
Net.fractionExcite = 0.5;
Net.ratioI = 3;
Net.R = 10;   % spectral radius

% create network
[W,Icells] = initialise_SOC_Weight_Matrix(Net.N,Net.p,Net.fractionExcite,Net.ratioI,Net.R);

% get eigenvalues of that network
initial_eigs = eig(W,'vector');
f = figure;
plot(real(initial_eigs),imag(initial_eigs),'.','Color',[0.7 0.7 0.7]); hold on

% save resulting network
fileID = datestr(now,30);

save(['Random_Network_' fileID],'W','Icells','Net');
