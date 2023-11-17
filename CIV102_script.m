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

ymaxs = zeros(1, numxsections);
for i = 1:numxsections
    ymaxs(i) = max(max(xsections{i}(:, 2)), max(xsections{i}(:, 4)));
end

glueheights = [75, 50];
numglues = 2;

%% 1. SFD, BMD under train loading
x_train = [52 228 392 568 732 908]; % Train Load Locations
x_train = x_train - 52;
P_train_LC1 = [1 1 1 1 1 1] * (P/6);
P_train_LC2 = [1.35 1.35 1 1 1 1] * (P/6);


P_train = P_train_LC2;

% train_locs = x(x < (L - x_train(end))) + 1;
train_locs = (-x_train(end)):1:L;

n_train = length(train_locs); % num of train locations
SFDi = zeros(n_train, n+1); % 1 SFD for each train loc.
BMDi = zeros(n_train, n+1); % 1 BMD for each train loc.

% Solve for SFD and BMD with the train at different locations
for i = 1:n_train
    % start location of train
    locs = train_locs(i) + x_train;

    logical = locs >= 0 & locs <= L;

    locs = locs(logical);

    
    p_applicable = P_train(logical);

    % sum of moments at A eqn
    By = sum((locs) .* p_applicable) / L;
    % sum of Fy eqn 
    Ay = sum(p_applicable) - By;
 
    % construct applied loads
    w = zeros(1, n + 1);
    w(locs + 1) = -p_applicable;
    w(1) = w(1) + Ay;
    w(L + 1) = w(L + 1) + By;
    
    SFDi(i, :) = cumsum(w);
    BMDi(i, :) = cumtrapz(SFDi(i,:));
    plot(x, SFDi(i, :), "r")
    hold on

end
hold off
SFE = max(abs(SFDi)); % SFD envelope
BME = max(BMDi); % BMD envelope

figure
plot(x, SFE)
figure
plot(x, BME)

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
        Iloc = b * h^3 / 12;
        ay2 = b * h * (abs(y - ybars(i)))^2;
        Itot = Itot + Iloc + ay2;
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



%% 3. Calculate Applied Stress
S_top = zeros(1, n + 1);
S_bot = zeros(1, n + 1);
T_cent = zeros(1, n + 1);
T_glue = zeros(numgles, n + 1);
for i = 1:n + 1
    S_top(i) = ymaxs(i) * BME(i) / I(i);
end
%% 4. Material and Thin Plate Buckling Capacities
E = 4000;
mu = 0.2;
S_tens = 
S_comp = 
T_max = 
T_gmax = 
S_buck1 = 
S_buck2 = 
S_buck3 = 
T_buck = 
%% 5. FOS
FOS_tens = 
FOS_comp = 
FOS_shear = 
FOS_glue = 
FOS_buck1 = 
FOS_buck2 = 
FOS_buck3 = 
FOS_buckV = 
%% 6. Min FOS and the failure load Pfail
minFOS = 
Pf =
%% 7. Vfail and Mfail
Mf_tens = 
Mf_comp = 
Vf_shear = 
Vf_glue = 
Mf_buck1 = 
Mf_buck2 = 
Mf_buck3 = 
Vf_buckV = 
%% 8. Output plots of Vfail and Mfail
subplot(2,3,1)
hold on; grid on; grid minor;
plot(x, Vf_shear, 'r')
plot(x, -Vf_shear.* SFD, 'r')
plot(x, SFDi, 'k');
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
legend('Matboard Shear Failure')
xlabel('Distance along bridge (mm)')
ylabel('Shear Force (N)')
