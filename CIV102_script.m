clear; close all;
%% 0. Initialize Parameters
L = 1200; % Length of bridge
n = 1200; % Discretize into 1 mm seg.
P = 400; % Total weight of train IF all are freight trains[N]
x = linspace(0, L, n+1); % x-axis generate n + 1 evenly spaced points

section1 = [100, 1.27, 1.27/2;
            80, 1.27, 1.27/2 + 75;
            1.27 * 2, 75 - 1.27, (75 - 1.27) / 2 + 1.27;
            10, 1.27, 3 * 1.27 / 2]; % b, h, ybar

section2 = [100, 1.27, 1.27/2;
            80, 1.27, 1.27/2 + 75;
            1.27 * 2, 75 - 1.27, (75 - 1.27) / 2 + 1.27;
            10, 1.27, 3 * 1.27 / 2]; % b, h, ybar

% section1 = [0, 75 + 1.27, 100, 75;
%             10, 75, 15 + 1.27, 75 - 1.27;
%             10, 75 - 1.27, 10 + 1.27, 0;
%             10 + 1.27, 1.27, 90 - 1.27, 0; % bottom center
%             85 - 1.27, 75, 90, 75 - 1.27;
%             90 - 1.27, 75 - 1.27, 90, 0];
% 
% section2 = [0, 75 + 1.27, 100, 75;
%             10, 75, 15 + 1.27, 75 - 1.27;
%             10, 75 - 1.27, 10 + 1.27, 0;
%             10 + 1.27, 1.27, 90 - 1.27, 0; % bottom center
%             85 - 1.27, 75, 90, 75 - 1.27;
%             90 - 1.27, 75 - 1.27, 90, 0];


xsections = {section1, section2}; % list of all cross sections

xsectionpts = [0, 1200]; % define where the x section change is. Size is n - 1

numxsections = length(xsections);



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
    w(round((locs + 1) / L * n)) = -p_applicable;
    w(1) = w(1) + Ay;
    w(round((L + 1.0) / L * n)) = w(round((L + 1.0) / L * n)) + By;
    
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

xsectionsdsc = cell(1, n + 1);
xsectionsdsc(round(xsectionpts(1) / L * n) + 1) = xsections(1); 
for i = 1:(numxsections - 1)
    idx1 = round(xsectionpts(i) / L * n) + 1;
    idx2 = round(xsectionpts(i + 1) / L * n) + 1;
    ratio = (xsections{i + 1} - xsections{i}) ./ (xsectionpts(i + 1) - xsectionpts(i));
    xsectionsdsc{idx1} = xsections{i};
    for j = (idx1 + 1):(idx2)
        xsectionsdsc{j} = xsections{i} + ratio * (j / n * L);
    end
end

ybots = zeros(1, n + 1);
for i = 1:n + 1
    ybots(i) = max(xsectionsdsc{i}(:, 2) / 2 + xsectionsdsc{i}(:, 3));
end

ybars = zeros(1, n + 1);
Is = zeros(1, n + 1);
QCent = zeros(1, n + 1);
Qglues = zeros(n + 1, numglues);



for i = 1:(n + 1)
    if i > 1 && isequal(xsectionsdsc{i}, xsectionsdsc{i - 1})
        ybars(i) = ybars(i - 1);
        Is(i) = Is(i - 1);
        QCent(i) = QCent(i - 1);
        Qglues(i) = Qglues(i - 1);
        continue
    end
    % ybar, I, Qcent
    % sum the areas times distance to centroid of each individual square
    sumAD = 0;
    sumAreas = 0;
    for j = 1:length(xsectionsdsc{i})
        % ybar
        b = xsectionsdsc{i}(j, 1);
        h = xsectionsdsc{i}(j, 2);
        y = xsectionsdsc{i}(j, 3);
        sumAD = sumAD + b * h * y;
        sumAreas = sumAreas + b * h;

    end
    ybars(i) = sumAD / sumAreas;

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
    Qgs = zeros(1, numglues);
    for j = 1:length(xsectionsdsc{i})
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
                Qc = Qc + area * (glueheights(j) - centY);
            end
        end
        Qgs(j) = Qtop;
    end
    QCent(i) = Qc;
    Qglues(i, :) = Qgs;
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
        y1 = y - h / 2;
        y2 = y + h / 2;
        if y2 >= ybars(i) && y1 < ybars(i)
            bBot = bBot + b;
        end
        if y2 > ybars(i) && y1 <= ybars(i)
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
            y1 = y - h / 2;
            y2 = y + h / 2;
            if y2 >= glueheights(j) && y1 < glueheights(j)
                bBot = bBot + b;
            end
            if y2 > glueheights(j) && y1 <= glueheights(j)
                bTop = bTop + b;
            end
        end
        bGluesTop(j, i) = bTop;
        bGluesBot(j, i) = bBot;
    end
end


%% 3. Calculate Applied Stress
S_top = zeros(1, n + 1);
S_bot = zeros(1, n + 1);
T_cent = zeros(1, n + 1);
T_glue = zeros(numglues, n + 1);
for i = 1:(n + 1)
    S_top(i) = -(ybots(i)- ybars(i)) * BME(i) / Is(i);
    S_bot(i) = ybars(i) * BME(i) / Is(i);
    T_cent(i) = SFE(i) * QCent(i) / Is(i) / min(bCentBot(i), bCentTop(i));
    for j = 1:numglues
        T_glue(j, i) = SFE(i) * QCent(i) / Is(i) / min(bGluesTop(j, i), bGluesBot(j, i));
    end
end

%% 4. Material and Thin Plate Buckling Capacities
% E = 4000;
% mu = 0.2;
% S_tens = 
% S_comp = 
% T_max = 
% T_gmax = 
% S_buck1 = 
% S_buck2 = 
% S_buck3 = 
% T_buck = 

%% 5. FOS
% FOS_tens =  
% FOS_comp = 
% FOS_shear = 
% FOS_glue = 
% FOS_buck1 = 
% FOS_buck2 = 
% FOS_buck3 = 
% FOS_buckV = 

%% 6. Min FOS and the failure load Pfail
% minFOS = 
% Pf =
% %% 7. Vfail and Mfail
% Mf_tens = 
% Mf_comp = 
% Vf_shear = 
% Vf_glue = 
% Mf_buck1 = 
% Mf_buck2 = 
% Mf_buck3 = 
% Vf_buckV = 


%% 8. Output plots of Vfail and Mfail
% subplot(2,3,1)
% hold on; grid on; grid minor;
% plot(x, Vf_shear, 'r')
% plot(x, -Vf_shear.* SFD, 'r')
% plot(x, SFDi, 'k');
% plot([0, L], [0, 0], 'k', 'LineWidth', 2)
% legend('Matboard Shear Failure')
% xlabel('Distance along bridge (mm)')
% ylabel('Shear Force (N)')
