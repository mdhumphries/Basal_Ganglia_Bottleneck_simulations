% script to test SC response to one input: where BG output and command signal
% coincide
%
% map network (x,y): no interactions between SC neurons
%
% Mark Humphries 

clearvars; close all;

saveID = 'GridOnly_SC_Dual_Input';
rng(1);  % set seed

%% SC grid parameters
Net.gMap = 20;  % number of nodes per side of (x,y) grid
Net.Nmap = Net.gMap ^2; % total number of grid nodes

% get every neuron's (x,y) co-ordinates
rep = repmat(1:Net.gMap,Net.gMap,1);
Net.ixCoords = [repmat([1:Net.gMap]',Net.gMap,1) rep(:)] ;  % array of all neuron co-ordinates

%% command input to grid

% local group of neurons with input
Sim.xyGroup1 = [5,9]; % set (x,y) location of group centre
Sim.dThresh = 2;  % distance threshold to be in group
Sim.Input1 = 10;  % input to group 1; 

% other group with input
Sim.xyGroup2 = [14, Sim.xyGroup1(2)]; % same column, different row
Sim.Input2 = 40;  % input to group 2; 


%% SNr parameters, including connections

% the number of SNrs
Net.SNr_xy_resolution = 1;  % number of coordinates per SNr neuron (minimum 1, maximum Net.gmap)
Net.SNr_per_Side = Net.gMap / Net.SNr_xy_resolution;   % on each side of the grid, up to a maxium of Net.gmap
Net.Nsnr = Net.SNr_per_Side * 2;  % total number of SNr neurons

% SNr input
Sim.bSNr = 30; % baseline activity: 30 spikes/s
Sim.tStep= 500; % duration of step changes (in time-steps)

% SNr neurons to turn off and/or on around (x,y) location of command input
Net.SNrRadius = floor(Sim.dThresh / Net.SNr_xy_resolution);   % radius of SNr neurons to turn off from centre point - counts along x or y axis
Sim.SNrOff = 0;  % what value is "off"
Sim.SNrOn = 40; % Sim.bSNr;  % what value is "on"  [40 to suppress second input; Sim.bSNr for no effect]

%% SC neuron parameters
% proportion of neurons receiving random input
Sim.Pinput = 0; % try 5 or 10% here
Sim.Minput = 5;
Sim.SDinput = 0;

% neuron properties
Sim.ts = 10; % ms
Sim.dt = 1; % ms
Sim.output = 'ReLu';

% duration of simulation
Sim.maxSteps = Sim.tStep * 3;
Sim.dRthresh = 0.1;  % stop once change in rate has fallen below this value...

%% SNr connections to SC: strip-tiling
Net.WSNr = Strip_tiling_2D(Net.gMap,Net.gMap,Net.SNr_per_Side);

%% superficial layer input connections to SC
Sim.Iext = zeros(Net.Nmap,1); 

% randomly selected neurons for background input
Ninput = round(Sim.Pinput * Net.Nmap);  % Number of SC neurons to give background input to
ids = randperm(Net.Nmap);  % randomised list of all neurons
Sim.Iext(ids(1:Ninput)) = Sim.SDinput*randn(Ninput,1) + Sim.Minput;  % assign that many neurons an input
Sim.Iext(Sim.Iext < 0) = 0;     % rectify

% group of neurons getting enhanced input
if ~isempty(Sim.xyGroup1)
    ix = sub2ind([Net.gMap Net.gMap],Sim.xyGroup1(1),Sim.xyGroup1(2));  % pick a neuron at location x,y
    dSrc = abs(Net.ixCoords(:,1) - Net.ixCoords(ix,1)) + abs(Net.ixCoords(:,2) - Net.ixCoords(ix,2));
    ids = find(dSrc <= Sim.dThresh);
    Sim.Iext(sub2ind([Net.gMap Net.gMap],Net.ixCoords(ids,1),Net.ixCoords(ids,2))) = Sim.Input1;
end 

% second group of neurons getting enhanced input
if ~isempty(Sim.xyGroup2)
    ix = sub2ind([Net.gMap Net.gMap],Sim.xyGroup2(1),Sim.xyGroup2(2));  % pick a neuron at location x,y
    dSrc = abs(Net.ixCoords(:,1) - Net.ixCoords(ix,1)) + abs(Net.ixCoords(:,2) - Net.ixCoords(ix,2));
    ids = find(dSrc <= Sim.dThresh);
    Sim.Iext(sub2ind([Net.gMap Net.gMap],Net.ixCoords(ids,1),Net.ixCoords(ids,2))) = Sim.Input2;
end 

%% SNr output

Sim.rSNr = zeros(Net.Nsnr,Sim.maxSteps) + Sim.bSNr; % all SNr neurons have baseline output

% construct OFF input 
% convert command input centre (x,y) to nearest SNr neuron
x_SNr = ceil(Sim.xyGroup1(1) / Net.SNr_xy_resolution);
y_SNr = Net.SNr_per_Side + ceil(Sim.xyGroup1(2) / Net.SNr_xy_resolution); 

% create set of all SNr neurons to turn off
x_SNr = x_SNr - Net.SNrRadius:x_SNr + Net.SNrRadius;
y_SNr = y_SNr - Net.SNrRadius:y_SNr + Net.SNrRadius;

% assign these neurons to be off for the duration of tStep
Sim.rSNr(x_SNr,Sim.tStep:2*Sim.tStep-1) = Sim.SNrOff;  % turn off one location
Sim.rSNr(y_SNr,Sim.tStep:2*Sim.tStep-1) = Sim.SNrOff;  % turn off one location

% construct ON input
if ~isempty(Sim.xyGroup2)
    x_SNr = ceil(Sim.xyGroup2(1) / Net.SNr_xy_resolution);  % x location of second input
    x_SNr = x_SNr - Net.SNrRadius:x_SNr + Net.SNrRadius;
    % assign these SNr neurons to be enhanced for the duration of tStep
    Sim.rSNr(x_SNr,Sim.tStep:2*Sim.tStep-1) = Sim.SNrOn;  % turn up one location
end

%% run model
Out.a = zeros(Net.Nmap,Sim.maxSteps);
Out.r = zeros(Net.Nmap,Sim.maxSteps);
I = zeros(Net.Nmap,1);
for t = 2:Sim.maxSteps
    % input
    I = Net.WSNr * Sim.rSNr(:,t-1) + Sim.Iext;
    
    % update activity
    Out.a(:,t) = Out.a(:,t-1) + Sim.dt*(-Out.a(:,t-1) + I)/Sim.ts;
    
    % get output
    Out.r(:,t) = neuron_output(Out.a(:,t),Sim.output);
    
    Out.dr(t) = sum(abs(Out.r(:,t) - Out.r(:,t-1)));
%     if Out.dr(t) < Sim.dRthresh & t > Sim.tStep+50  % if settled after change to input, then quit
%         
%         break
%     end
end
Out.nSteps = t;

%% save parameters, network, and output
if exist('saveID','var')
    save(saveID,'Net','Sim','Out');
end

%% plot final time-point on grid

% reshape is column by column
matI = reshape(Sim.Iext,Net.gMap,Net.gMap);

% take snapshots of SNr inputs and SC outputs
snappts = (Sim.tStep - 1):Sim.tStep:Sim.maxSteps; % take a snapshot of end point after every change
Nsnap = numel(snappts);
for iSnap = 1:Nsnap
    rblock = Out.r(:,round(snappts(iSnap)));
    rSnap(:,:,iSnap) = reshape(rblock,Net.gMap,Net.gMap);
    
    SNrinputBlock = Net.WSNr * Sim.rSNr(:,round(snappts(iSnap)));
    SNrInputSnap(:,:,iSnap) = reshape(SNrinputBlock,Net.gMap,Net.gMap); 
end

minrSnap = min(min(min(rSnap)));
maxrSnap = max(max(max(rSnap)));

minSNrSnap = min(min(min(SNrInputSnap)));
maxSNrSnap = max(max(max(SNrInputSnap)));


% external input
figure
imagesc(matI)
axis square
colorbar
title('External input')

% sequence of SNr inputs and SC outputs
figure
for iSnap = 1:Nsnap
    subplot(2,Nsnap,iSnap),imagesc(squeeze(SNrInputSnap(:,:,iSnap)));
    axis square
    ax = gca;
    ax.CLim = [minSNrSnap maxSNrSnap];
    colorbar
    title(['SNr input: ' num2str(round(snappts(iSnap)))])
    % plot these on same color scale
    subplot(2,Nsnap,Nsnap+iSnap),imagesc(squeeze(rSnap(:,:,iSnap)));
    axis square
    ax = gca;
    ax.CLim = [minrSnap maxrSnap];
    colorbar
    title(['SC output: ' num2str(round(snappts(iSnap)))])
end

% output
figure
plot(Out.r(:,1:Out.nSteps)')
xlabel('Time (ms)')
ylabel('Rate (spikes/s)')





