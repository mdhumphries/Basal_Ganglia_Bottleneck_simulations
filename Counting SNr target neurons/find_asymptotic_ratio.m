function [data_asymptote,fit_asymptote,datafit] = find_asymptotic_ratio(x,y,eta)

% FIND_ASYMPTOTIC_RATIO finds the asymptotic value of data
% [DA,FA,FIT] = FIND_ASYMPTOTIC_RATIO(X,Y,ETA) finds the asymptotic value of the
% data values Y that are a function of X. The threshold ETA defines the
% size of step change between adjacent values of Y that are considered
% asymptotic
%
% DA: asymptote estimated directly from data 
% FA: asymptote estimated from polynomial fit to data
% Both check for two consecutive data differences that are less than ETA
% Values will be NaN no threshold crossing is found
% FIT: the fitted model (a FIT object - see "fit" for details)
%
% Mark Humphries 10/5/2024

if numel(x) ~= numel(y) 
    error('x and y must have the same number of entries');
end

% make column vectors
if size(x,1) == 1 x = x'; end
if size(y,1) == 1 y = y'; end

% defaults
data_asymptote = nan;
fit_asymptote = nan;

% option 1: direct data estimate
dy = diff(y);
crossing = dy <= eta;
ix_consecutive = find(diff(crossing) == 0,1,'last');
if ~isempty(ix_consecutive)  
    data_asymptote = y(ix_consecutive+2);
end

% option 2: fit model
datafit = fit(x,y,'exp2');
dx = mean(diff(x));  % step size of data
x_eval = (min(x)-dx):dx*0.1:(max(x)+dx);  % data points to evaluate fitted curve
dy = differentiate(datafit,x_eval);       % differentiate curve - has value for every point on x  
crossing = dy <= eta;                     % find its asymptote
ix_consecutive = find(diff(crossing) == 0,1,'last'); % check consecutive points are below threshold
if ~isempty(ix_consecutive) 
    yfit = feval(datafit,x_eval);         % evaluate fitted curve
    try
    fit_asymptote = yfit(ix_consecutive+1); % find asymptotic value: 
    catch
        keyboard
    end
end
figure
plot(datafit,x,y)

