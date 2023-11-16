clear; close all;
%% 0. Initialize Parameters
L = 1200; % Length of bridge
n = 1200; % Discretize into 1 mm seg.
P = 400; % Total weight of train IF all are freight trains[N]
x = linspace(0, L, n+1); % x-axis generate n + 1 evenly spaced points

xsection1 = [0, 75 + 1.27, 100, 75;
            10, 75, 15 + 1.27, 75 - 1.27;
            10, 75 - 1.27, 10 + 1.27, 0;
            10 + 1.27, 1.27, 90 - 1.27, 0; % bottom center
            85 - 1.27, 75, 90, 75 - 1.27;
            90 - 1.27, 75 - 1.27, 90, 0];
xsection2 = [0, 75 + 1.27, 100, 75;
            10, 75, 15 + 1.27, 75 - 1.27;
            10, 75 - 1.27, 10 + 1.27, 0;
            10 + 1.27, 1.27, 90 - 1.27, 0; % bottom center
            85 - 1.27, 75, 90, 75 - 1.27;
            90 - 1.27, 75 - 1.27, 90, 0];

xsections = {xsection1, xsection2}; % list of all cross sections

xsectionpts = [0, 1200]; % define where the x section change is. Size is n - 1

numxsections = length(xsections);

glueheights = [75, 50];
numglues = 2;

%% 1. SFD, BMD under train loading
x_train = [52 228 392 568 732 908] - 51; % Train Load Locations
P_train_LC1 = [1 1 1 1 1 1] * P/6;
P_train_LC2 = [1.35 1.35 1 1 1 1] * P;

P_train = P_train_LC2;

train_locs = x(x < (L - x_train(end))) + 1;
train_locs = (1-x_train(end)):1:L - 2;

n_train = length(train_locs); % num of train locations
SFDi = zeros(n_train, n+1); % 1 SFD for each train loc.
BMDi = zeros(n_train, n+1); % 1 BMD for each train loc.

% Solve for SFD and BMD with the train at different locations
for i = 1:n_train
    % start location of train
    locs = train_locs(i) + x_train;

    locs = locs(locs > 1 & locs < L - 2);
    
    p_applicable = P_train(locs > 1 & locs < L - 2);

    % sum of moments at A eqn
    By = sum((locs) .* p_applicable) / L;
    % sum of Fy eqn 
    Ay = sum(p_applicable) - By;
 
    % construct applied loads
    w = zeros(1, n + 1);
    w(locs) = -p_applicable;
    w(1) = Ay;
    w(n + 1) = By;
    
    SFDi(i, :) = cumsum(w);
    BMDi(i, :) = cumsum(SFDi(i,:));
    plot(x + 1, SFDi(i, :))
    hold on

end
hold off
SFE = max(abs(SFDi)); % SFD envelope
BME = max(BMDi); % BMD envelope

figure
plot(x + 1, SFE)
figure
plot(x + 1, BME)

%% 2. Geometric properties of cross sections



ybars = zeros(1, numxsections);
Is = zeros(1, numxsections);
QCent = zeros(1, numxsections);
Qglues = zeros(numxsections, numglues);

for i = 1:numxsections
    % ybar, I, Qcent
    % sum the areas times distance to centroid of each individual square
    sumAD = 0;
    sumAreas = 0;
    for j = 1:length(xsections{i})
        % ybar
        ax = xsections{i}(j, 1);
        ay = xsections{i}(j, 2);
        bx = xsections{i}(j, 3);
        by = xsections{i}(j, 4);
        sumAD = sumAD + (bx - ax) * (ay - by) * (by + ay) / 2;
        sumAreas = sumAreas + (bx - ax) * (ay - by);

    end
    ybars(i) = sumAD / sumAreas;

    % I
    % bh^3/12 for every bit
    Itot = 0;
    for j = 1:length(xsections{i})
        ax = xsections{i}(j, 1);
        ay = xsections{i}(j, 2);
        bx = xsections{i}(j, 3);
        by = xsections{i}(j, 4);
        y = (by + ay) / 2;
        h = (ay - by);
        b = (bx - ax);
        I = b * h^3 / 12;
        ay2 = b * h * (abs(y - ybars(i)))^2;
        Itot = Itot + I + ay2;
    end
    Is(i) = Itot;
    
    % all Qs
    Qc = 0;
    Qgs = zeros(1, numglues);
    for j = 1:length(xsections{i})
        % calculate area above centroidal axis of the jth shape, as well as
        % centroid position relative to y axis
        ax = xsections{i}(j, 1);
        ay = xsections{i}(j, 2);
        bx = xsections{i}(j, 3);
        by = xsections{i}(j, 4);
        if ay > ybars(i)
            area = (bx - ax) * (ay - ybars(i));
            centY = (ay + ybars(i)) / 2;
            Qc = Qc + area * (centY - ybars(i));
        end
    end
    for j = 1:numglues
        % calculate Q from top
        Qtop = 0;
        for k = 1:length(xsections{i})
            ax = xsections{i}(j, 1);
            ay = xsections{i}(j, 2);
            bx = xsections{i}(j, 3);
            by = xsections{i}(j, 4);
            if ay > glueheights(j)
                area = (bx - ax) * (ay - glueheights(j));
                centY = (ay + glueheights(j)) / 2;
                Qtop = Qtop + area * (centY - glueheights(j));
            end
        end
        Qgs(j) = Qtop;
    end
    QCent(i) = Qc;
    Qglues(i, :) = Qgs;
end


