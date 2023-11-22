function [fails, totLoad, FOSs] = checkPFail(P)
    %% 0. Initialize Parameters
    L = 1200; % Length of bridge
    n = 1200; % Discretize into 1 mm seg.
    x = linspace(0, L, n+1); % x-axis generate n + 1 evenly spaced points
    
    % Cross sections are defined by an nx3 matrix. Each row represents a
    % rectangle. 
    % section[i,1] = b, section[i, 2] = h, section[i, 3] = ybar
    % **Coordinates measured as right and down positive**
    section1 = [100, 1.27, 1.27/2;
                80, 1.27, 1.27/2 + 75;
                1.27 * 2, 75 - 1.27, (75 - 1.27) / 2 + 1.27;
                10, 1.27, 3 * 1.27 / 2]; % b, h, ybar
    
    xsections = {section1, section1}; % cell array of all cross sections
    
    % define the location of each cross section. At least two required.
    xsectionpts = [0, 1200]; 
    numxsections = length(xsections);
    
    % Define top flange parameters
    topConstThick = 1.27;
    topConstWidth = 80 - 2 * 1.27;
    topFreeThick = 1.27;
    topFreeWidth = 10;
    
    % Define glue joints
    glueheights = [1.27]; %#ok<*NBRAK2>
    numglues = size(glueheights, 1);
    
    %% 1. SFD, BMD under train loading
    x_train = [52 228 392 568 732 908]; % Train Load Locations
    x_train = x_train - 52;
    P_train_LC1 = [1 1 1 1 1 1] * (P/6);
    P_train_LC2 = [1.35 1.35 1 1 1 1] * (P/6);
    
    % **Modify to change load case**
    P_train = P_train_LC2;

    totLoad = sum(P_train);
    
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
    hold off
    SFE = max(abs(SFDi)); % SFD envelope
    BME = max(BMDi); % BMD envelope
    
    %% 2. Geometric properties of cross sections
    
    % Discretized cross-sections. One cross section per mm.
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
    Areas = zeros(1, n + 1);
    
    % ybar, I, Qcent, Qglue
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
    T_glue = zeros(numglues, n + 1);
    S_top = -ybars .* BME ./ Is;
    S_bot = (ybots - ybars) .* BME ./ Is;
    T_cent = SFE .* QCent ./ Is ./ min(bCentBot, bCentTop);
    
    for j = 1:numglues
        T_glue(j, :) = SFE .* QCent ./ Is ./ min(bGluesTop(j, :), bGluesBot(j, :));
    end
    
    %% 4. Material and Thin Plate Buckling Capacities
    diaphragmDist = 1200; % distance between diaphragms. Assume constant for ease of calculation.
    webHs = zeros(1, n + 1);
    for i = 1:n + 1
        cnt = 1;
        for j = 1:size(xsectionsdsc{i}, 1)
            b = xsectionsdsc{i}(j, 1);
            h = xsectionsdsc{i}(j, 2);
            y = xsectionsdsc{i}(j, 3);
            if h > b
                webHs(cnt) = max(webHs(cnt), xsectionsdsc{i}(j, 2));
                cnt = cnt + 1;
            end
        end
    end
    
    E = 4000;
    mu = 0.2;
    t = 1.27;
    S_tens = max(S_top, S_bot);
    S_comp = abs(min(S_top, S_bot));
    T_max = T_cent;
    T_gmax = max(T_glue);
    S_buck1 = 4*pi^2*E ./ (12 * (1 - mu^2)) .* ((topConstThick/topConstWidth).^2);
    S_buck2 = 0.4254*pi^2*E ./ (12 * (1 - mu^2)) .* ((topFreeThick/topFreeWidth).^2);
    
    S_buck3 = 6*pi^2*E ./ (12*(1-mu^2)) .* (t./(ybars - topConstThick)).^2;
    
    T_buck = 5*pi^2*E ./ (12*(1-mu^2)) .* ((t./webHs).^2 + (t/diaphragmDist).^2);
    
    
    %% 5. FOS
    FOS_tens = 30 ./ S_tens;
    FOS_comp = 6 ./ S_comp;
    FOS_shear = 4 ./ T_max;
    FOS_glue = 2 ./ T_gmax;
    FOS_buck1 = S_buck1 ./ S_comp;
    FOS_buck2 = S_buck2 ./ S_comp;
    FOS_buck3 = S_buck3 ./ S_comp;
    FOS_buckV = T_buck ./ (SFE / Areas);

    FOSs = [min(FOS_tens), min(FOS_comp), min(FOS_shear), min(FOS_glue), min(FOS_buck1), min(FOS_buck2), min(FOS_buck3), min(FOS_buckV)];
    
    %% 6. Min FOS and the failure load Pfail
    minFOS = min(FOSs);

    if minFOS < 1
        fails = true;
    else
        fails = false;
    end
end