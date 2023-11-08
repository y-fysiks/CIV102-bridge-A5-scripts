clear; close all;
%% 0. Initialize Parameters
L = 1200; % Length of bridge
n = 1200; % Discretize into 1 mm seg.
P = 400; % Total weight of train [N]
x = linspace(0, L, n+1); % x-axis generate n + 1 evenly spaced points

xsection1 = [0, 75 + 1.27, 100, 75;
            10, 75, 15 + 1.27, 75 - 1.27;
            10, 75 - 1.27, 10 + 1.27, 0;
            10 + 1.27, 1.27, 90 - 1.27, 0; % bottom center
            85 - 1.27, 75, 90, 75 - 1.27;
            90 - 1.27, 75 - 1.27, 90, 0];

xsections = [xsection1]; % list of all cross sections

xsectionchanges = []; % define where the x section change is. Size is n - 1

numxsections = length(xsections);

%% 1. SFD, BMD under train loading
x_train = [52 228 392 568 732 908]; % Train Load Locations
P_train = [1 1 1 1 1 1] * P/6;

train_locs = [480, 0];
n_train = length(train_locs); % num of train locations
SFDi = zeros(n_train, n+1); % 1 SFD for each train loc.
BMDi = zeros(n_train, n+1); % 1 BMD for each train loc.

% Solve for SFD and BMD with the train at different locations
for i = 1:n_train 
    % start location of train
    loc = train_locs(i);
    
    % sum of moments at A eqn
    % sum of Fy eqn 
 
    % construct applied loads
    % w(x)
 
    % SFD = num. integral(w)
    % BMD = num. integral(SFD)
end
SFD = max(abs(SFDi)); % SFD envelope
BMD = max(BMDi); % BMD envelope