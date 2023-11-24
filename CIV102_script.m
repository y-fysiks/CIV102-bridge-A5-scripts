clear; close all;
%% 0. Initialize Parameters
L = 1200; % Length of bridge
n = 1200; % Discretize into 1 mm seg.
P = 400; % Total weight of train IF all are freight trains[N]
x = linspace(0, L, n+1); % x-axis generate n + 1 evenly spaced points

% Cross sections are defined by an nx3 matrix. Each row represents a
% rectangle. 
% section[i,1] = b, section[i, 2] = h, section[i, 3] = ybar
% **Coordinates measured as right and down positive**

%% design 0
section1 = [100, 1.27, 1.27/2;
            80, 1.27, 1.27/2 + 75;
            1.27 * 2, 75 - 1.27, (75 - 1.27) / 2 + 1.27;
            10, 1.27, 3 * 1.27 / 2;
            10, 0.0001, 1.27]; % b, h, ybar. last line is glue width,

% Iteration 1
sectionv2 = [100, 1.27*2, 1.27; % top flange
            80, 1.27, 1.27/2 + 75 + 1.27; % bottom flange
            1.27 * 2, 75 - 1.27, (75 - 1.27) / 2 + 2*1.27; % 
            10, 1.27, 2*1.27 + 1.27 / 2; % glue tabs
            10, 0.0001, 2*1.27]; % b, h, ybar. last line is glue width,

% Iteration 2
sectionv3 = [100, 1.27*2, 1.27; % top flange
            1.27 * 2, 75, (75) / 2 + 2*1.27; % 
            10, 1.27, 2*1.27 + 1.27 / 2; % glue tabs
            10, 0.0001, 2*1.27]; % b, h, ybar. last line is glue width,

% Iteration 3
sectionv4 = [110, 1.27*2, 1.27; % top flange
            1.27 * 2, 75, (75) / 2 + 2*1.27; % 
            10, 1.27, 2*1.27 + 1.27 / 2; % glue tabs
            10, 0.0001, 2*1.27]; % b, h, ybar. last line is glue width,

% Iteration 4
section1v5 = [110, 1.27*2, 1.27; % top flange
            1.27 * 2, 75, 75 / 2 + 2*1.27; % 
            10, 1.27, 2*1.27 + 1.27 / 2; % glue tabs
            10, 0.0001, 2*1.27]; % b, h, ybar. last line is glue width,

section2v5 = [110, 1.27*2, 1.27; % top flange
            1.27 * 2, 150, 150 / 2 + 2*1.27; % 
            10, 1.27, 2*1.27 + 1.27 / 2; % glue tabs
            10, 0.0001, 2*1.27]; % b, h, ybar. last line is glue width,

% Iteration 5
section1v6 = [110, 1.27*2, 1.27; % top flange
            1.27 * 2, 120, 100 / 2 + 2*1.27; % 
            10, 1.27, 2*1.27 + 1.27 / 2; % glue tabs
            10, 0.0001, 2*1.27]; % b, h, ybar. last line is glue width,

section2v6 = [140, 1.27*2, 1.27; % top flange
            1.27 * 2, 160, 120 / 2 + 2*1.27; % 
            10, 1.27, 2*1.27 + 1.27 / 2; % glue tabs
            10, 0.0001, 2*1.27]; % b, h, ybar. last line is glue width,

% Iteration 6
section1v7 = [100, 1.27*2, 1.27; % top flange
            1.27 * 2, 120, 100 / 2 + 2*1.27; % 
            10, 1.27, 2*1.27 + 1.27 / 2; % glue tabs
            10, 0.0001, 2*1.27]; % b, h, ybar. last line is glue width,

section2v7 = [100, 1.27*2, 1.27; % top flange
            1.27 * 2, 160, 120 / 2 + 2*1.27; % 
            10, 1.27, 2*1.27 + 1.27 / 2; % glue tabs
            10, 0.0001, 2*1.27]; % b, h, ybar. last line is glue width,

xsections = {section1v7, section2v7, section2v7, section1v7}; % cell array of all cross sections


% define the location of each cross section. At least two required.
xsectionpts = [0, 500, 700, 1200]; 

% Define top flange parameters
topConstThick = 1.27*2;
topConstWidth = 80;
topFreeThick = 1.27*2;
topFreeWidth = 10;

% Define glue joints
glueheights = [2*1.27, 1.27];

% distance between diaphragms. Assume constant for ease of calculation.
diaphragmDist = 100;

numglues = length(glueheights);

numxsections = length(xsections);

%% 1. SFD, BMD under train loading
x_train = [52 228 392 568 732 908]; % Train Load Locations
x_train = x_train - 52;
P_train_LC1 = [1 1 1 1 1 1] * (P/6);
P_train_LC2 = [1.35 1.35 1 1 1 1] * (P/6);

P_train = P_train_LC1;

% array of locations the end of the train can take, in intervals of 1mm
train_locs = (-x_train(end)):1:L;

n_train = length(train_locs); % number of train locations
SFDi = zeros(n_train, n+1); % 1 SFD for each train location
BMDi = zeros(n_train, n+1); % 1 BMD for each train location

% Solve for SFD and BMD with the train at different locations
for i = 1:n_train
    % start location of train
    locs = train_locs(i) + x_train;

    % ignore locations not on bridge
    logical = locs >= 0 & locs <= L; 
    locs = locs(logical);
    p_applicable = P_train(logical); % ignore corresponding loads

    % sum of moments at A eqn
    By = sum((locs) .* p_applicable) / L;

    % sum of Fy eqn 
    Ay = sum(p_applicable) - By;
 
    % construct applied loads
    w = zeros(1, n + 1);
    w(round((locs + 1) / L * n)) = -p_applicable;
    w(1) = w(1) + Ay;
    w(round((L + 1.0) / L * n)) = w(round((L + 1.0) / L * n)) + By;
    
    % Integrate the applied loads to get SFD
    SFDi(i, :) = cumsum(w);
    % Integrate the SFD to get BMD
    BMDi(i, :) = cumsum(SFDi(i,:));

end
SFE = max(abs(SFDi)); % SFD envelope
BME = max(BMDi); % BMD envelope

figure
subplot(1, 2, 1)
plot(x, SFE);


%% 2. Geometric properties of cross sections

% Discretized cross-sections. One cross section per mm.
xsectionsdsc = cell(1, n + 1);
xsectionsdsc(round(xsectionpts(1) / L * n) + 1) = xsections(1); 
for i = 1:(numxsections - 1)
    idx1 = round(xsectionpts(i) / L * n) + 1;
    idx2 = round(xsectionpts(i + 1) / L * n) + 1;
    ratio = (xsections{i + 1} - xsections{i}) ./ (idx2 - idx1);
    xsectionsdsc{idx1} = xsections{i};
    xsectionsdsc{idx2} = xsections{i + 1};
    for j = 1:(idx2-idx1-1)
        xsectionsdsc{idx1 + j} = xsections{i} + ratio * (j);
    end
end

for i = 1:n + 1
    
end

ybots = zeros(1, n + 1);
ybars = zeros(1, n + 1);
Is = zeros(1, n + 1);
QCent = zeros(1, n + 1);
Qglues = zeros(numglues, n + 1);
Areas = zeros(1, n + 1);

% ybar, I, Qcent, Qglues
for i = 1:(n + 1)
    % if i > 1 && isequal(xsectionsdsc{i}, xsectionsdsc{i - 1})
    %     ybars(i) = ybars(i - 1);
    %     Is(i) = Is(i - 1);
    %     QCent(i) = QCent(i - 1);
    %     Qglues(:, i) = Qglues(:, i - 1);
    %     Areas(i) = Areas(i - 1);
    %     continue
    % end
    % ybar, I, Qcent
    % sum the areas times distance to centroid of each individual square
    sumAD = 0;
    sumAreas = 0;
    for j = 1:size(xsectionsdsc{i}, 1)
        % ybar
        b = xsectionsdsc{i}(j, 1);
        h = xsectionsdsc{i}(j, 2);
        y = xsectionsdsc{i}(j, 3);
        sumAD = sumAD + b * h * y;
        sumAreas = sumAreas + b * h;
        Areas(i) = Areas(i) + b * h;
    end
    ybars(i) = sumAD / sumAreas;

    ybots(i) = max(xsectionsdsc{i}(:, 2) ./ 2 + xsectionsdsc{i}(:, 3));

    % I
    % bh^3/12 for every bit
    Itot = 0;
    for j = 1:size(xsectionsdsc{i}, 1)
        b = xsectionsdsc{i}(j, 1);
        h = xsectionsdsc{i}(j, 2);
        y = xsectionsdsc{i}(j, 3);
        Iloc = b * h^3 / 12;
        ay2 = b * h * (abs(y - ybars(i)))^2;
        Itot = Itot + Iloc + ay2;
    end
    Is(i) = Itot;
    
    % all Qs
    Qc = 0;
    for j = 1:size(xsectionsdsc{i}, 1)
        % calculate area above centroidal axis of the jth shape, as well as
        % centroid position relative to y axis
        b = xsectionsdsc{i}(j, 1);
        h = xsectionsdsc{i}(j, 2);
        y = xsectionsdsc{i}(j, 3);
        y1 = y - h/2;
        y2 = y + h/2;
        if y1 < ybars(i)
            area = b * (min(y2, ybars(i)) - y1);
            centY = (min(y2, ybars(i)) + y1) / 2;
            Qc = Qc + area * (ybars(i) - centY);
        end
    end
    for j = 1:numglues
        % calculate Q from top
        Qtop = 0;
        for k = 1:size(xsectionsdsc{i}, 1)
            % calculate area above centroidal axis of the jth shape, as well as
            % centroid position relative to y axis
            b = xsectionsdsc{i}(k, 1);
            h = xsectionsdsc{i}(k, 2);
            y = xsectionsdsc{i}(k, 3);
            y1 = y - h/2;
            y2 = y + h/2;
            if y1 < glueheights(j)
                area = b * (min(y2, glueheights(j)) - y1);
                centY = (min(y2, glueheights(j)) + y1) / 2;
                Qtop = Qtop + area * (ybars(i) - centY);
            end
        end
        Qglues(j, i) = Qtop;
    end
    QCent(i) = Qc;
end

% calculate max b
bCentTop = zeros(1, n + 1);
bCentBot = zeros(1, n + 1);

bGluesTop = zeros(numglues, n + 1);
bGluesBot = zeros(numglues, n + 1);


for i = 1:n + 1
    bBot = 0;
    bTop = 0;
    for j = 1:size(xsectionsdsc{i}, 1)
        b = xsectionsdsc{i}(j, 1);
        h = xsectionsdsc{i}(j, 2);
        y = xsectionsdsc{i}(j, 3);
        yb = round(ybars(i), 2);
        y1 = round(y - (h / 2), 2);
            y2 = round(y + (h / 2), 2);
        if y2 >= yb && y1 < yb
            bBot = bBot + b;
        end
        if y2 > yb && y1 <= yb
            bTop = bTop + b;
        end
    end
    bCentTop(i) = bTop;
    bCentBot(i) = bBot;
    for j = 1:numglues
        bBot = 0;
        bTop = 0;
        for k = 1:size(xsectionsdsc{i}, 1)
            b = xsectionsdsc{i}(k, 1);
            h = xsectionsdsc{i}(k, 2);
            y = xsectionsdsc{i}(k, 3);
            y1 = round(y - (h / 2), 5);
            y2 = round(y + (h / 2), 5);
            glueh = round(glueheights(j), 5);
            if y2 > glueh && y1 < glueh
                bBot = bBot + b;
            end

            if y2 > glueh && y1 <= glueh
                bTop = bTop + b;
            end
        end
        bGluesTop(j, i) = bTop;
        bGluesBot(j, i) = bBot;
    end
end


%% 3. Calculate Applied Stress
T_glue = zeros(numglues, n + 1);
S_top = ybars .* BME ./ Is;
S_bot = (ybots - ybars) .* BME ./ Is;
T_cent = SFE .* QCent ./ Is ./ min(bCentBot, bCentTop);

for j = 1:numglues
    T_glue(j, :) = SFE .* Qglues(j, :) ./ Is ./ min(bGluesTop(j, :), bGluesBot(j, :));
end

%% 4. Material and Thin Plate Buckling Capacities

E = 4000;
mu = 0.2;
t = 1.27;

S_tens = S_bot;
S_comp = S_top;
T_max = T_cent;
T_gmax = max(T_glue);
% Buckling case 1 (middle of top flange)
S_buck1 = 4*pi^2*E ./ (12 * (1 - mu^2)) .* ((topConstThick/topConstWidth)^2);
% Buckling case 2 (sides of top flange)
S_buck2 = 0.4254*pi^2*E ./ (12 * (1 - mu^2)) .* ((topFreeThick/topFreeWidth)^2);
% Buckling case 3 (Webs)
S_buck3 = 6*pi^2*E ./ (12*(1-mu^2)) .* (t./(ybars - topConstThick)).^2;
% Shear buckling
T_buck = 5*pi^2*E ./ (12*(1-mu^2)) .* (((t./(ybots - topConstThick)).^2) + ((t/diaphragmDist).^2));


%% 5. FOS
FOS_tens = 30 ./ S_tens;
FOS_comp = 6 ./ S_comp;
FOS_shear = 4 ./ T_max;
FOS_glue = 2 ./ T_gmax;
FOS_buck1 = S_buck1 ./ S_comp;
FOS_buck2 = S_buck2 ./ S_comp;
FOS_buck3 = S_buck3 ./ S_comp;
FOS_buckV = T_buck ./ abs(T_max);

%% 6. Min FOS and the failure load Pfail
FOSs = [min(FOS_tens), min(FOS_comp), min(FOS_shear), min(FOS_glue), min(FOS_buck1), min(FOS_buck2), min(FOS_buck3), min(FOS_buckV)];
minFOS = min(FOSs);
Pf = sum(P_train) * minFOS;

fprintf(['Factors of Safety\n' ...
         '  Matboard Tensile: %.2f\n' ...
         '  Matboard Compressive: %.2f\n' ...
         '  Matboard Shear: %.2f\n' ...
         '  Glue Shear: %.2f\n' ...
         '  Buckling Case 1: %.2f\n' ...
         '  Buckling Case 2: %.2f\n' ...
         '  Buckling Case 3: %.2f\n' ...
         '  Shear Buckling: %.2f\n'], FOSs)
fprintf('\nFailure load: %.2f N\n', Pf)
%% 7. Vfail and Mfail
Mf_tens = FOS_tens .* BME;
Mf_comp = FOS_comp .* BME;
Vf_shear = FOS_shear .* SFE;
Vf_glue = FOS_glue .* SFE;
Mf_buck1 = FOS_buck1 .* BME;
Mf_buck2 = FOS_buck2 .* BME;
Mf_buck3 = FOS_buck3 .* BME;
Vf_buckV = FOS_buckV .* SFE;


%% 8. Output plots of Vfail and Mfail
subplot(2,3,1)
hold on; grid on; grid minor;
plot(x, Vf_shear, 'r')
plot(x, SFE, 'k');
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
legend('Matboard Shear Failure')
xlabel('Distance along bridge (mm)')
ylabel('Absolute Max. Shear Force (N)')

subplot(2, 3, 2)
hold on; grid on; grid minor;
plot(x, Vf_glue, 'r')
plot(x, SFE, 'k');
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
legend('Glue Shear Failure')
xlabel('Distance along bridge (mm)')
ylabel('Absolute Max. Shear Force (N)')

subplot(2, 3, 3)
hold on; grid on; grid minor;
plot(x, Vf_buckV, 'r')
plot(x, SFE, 'k');
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
legend('Matboard Shear Buckling Failure')
xlabel('Distance along bridge (mm)')
ylabel('Absolute Max. Shear Force (N)')

subplot(2, 3, 4)
hold on; grid on; grid minor;
plot(x, Mf_tens, 'r')
plot(x, Mf_comp, 'b')
plot(x, BME, 'k');
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
legend('Matboard Tension Failure', 'Matboard Compression Failure')
xlabel('Distance along bridge (mm)')
ylabel('Bending moment (Nmm)')

subplot(2, 3, 5)
hold on; grid on; grid minor;
plot(x, Mf_buck1, 'r')
plot(x, Mf_buck2, 'b')
plot(x, BME, 'k');
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
legend('Matboard Buckling Failure, Top Flange - Mid', 'Matboard Buckling Failure, Top Flange - Sides')
xlabel('Distance along bridge (mm)')
ylabel('Bending moment (Nmm)')

subplot(2, 3, 6)
hold on; grid on; grid minor;
plot(x, Mf_buck3, 'r')
plot(x, BME, 'k');
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
legend('Matboard Buckling Failure, Webs')
xlabel('Distance along bridge (mm)')
ylabel('Bending moment (Nmm)')
