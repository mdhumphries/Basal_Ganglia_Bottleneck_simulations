function W = Strip_tiling_2D(nRows,nColumns,nInputs_per_side)

% STRIP_TILING_2D create strip-tiled basis functions of 2D space
% W = STRIP_TILING_2D(ROWS,COLS,INPUTS_PER_SIDE) creates an adjacency 
% matrix that defines the row-by-column strip-tiling of inputs to a grid. 
% The grid is ROWS x COLS nodes, and each axis receives INPUTS_PER_SIDE, 
% giving 2 x INPUTS_PER_SIDE inputs in total.
%
% W is an adjacency matrix of size [(ROWSxCOLS), 2xINPUTS_PER_SIDE], each
% row of W giving the inputs to that node. The inputs are strip-tiled: one
% input per row, one input per column.
%
% Feb 2016: Initial code
% 10/3/2022: created function, simplified code
% Mark Humphries 

total_number_of_target_coordinates = nRows * nColumns;
total_inputs = 2 * nInputs_per_side;
W = zeros(total_number_of_target_coordinates,total_inputs);  % to,from

% strip tile rows
rowwidth = nRows / nInputs_per_side;  % width of 1 strip, in rows
dcolstrip = 1:nColumns; % indices is all columns
strpctr = 1;
for iS = 1:nInputs_per_side
    indices_of_rowstrip = (strpctr-1)*rowwidth+1:strpctr*rowwidth; % indices are next strip of rows
    targetIDs = (repmat(dcolstrip,rowwidth,1)-1).*nRows + repmat(indices_of_rowstrip',1,nColumns);  % combine row and column sets to get complete indices
    W(targetIDs,iS) = -1;
    % keyboard
    strpctr = strpctr + 1;
end

% strip tile columns       
colwidth = nColumns/nInputs_per_side;
indices_of_rowstrip = 1:nRows;
strpctr = 1;

for iS = nInputs_per_side+1:2*nInputs_per_side    
    dcolstrip = (strpctr-1)*colwidth+1:strpctr*colwidth;
    targetIDs = (repmat(dcolstrip,nRows,1)-1).*nRows + repmat(indices_of_rowstrip',1,colwidth);
    W(targetIDs,iS) = -1;
    strpctr = strpctr + 1;   
end

%figure
%imagesc(W)

