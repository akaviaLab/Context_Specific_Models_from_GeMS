function runThresholdsINIT(figName, bb, modelName, cellLine, overWrite)

if (nargin < 5 || isempty(overWrite))
    overWrite = true;
end

%initCobraToolbox
load(['ID_FPKM_', cellLine, '.mat'], 'num');
load(['growthRate_',cellLine,'.mat'], 'blb')
load(['gene_threshold_',cellLine,'.mat'], 'ths')

if (isempty(modelName))
    load(['model_u_',cellLine,'.mat'], 'model_u')
    load(['model_c_',cellLine,'.mat'], 'model_c')
    load(['model_s_',cellLine,'.mat'], 'model_s')
else
    load([modelName, '_', cellLine, '.mat'], 'model_u','model_c', 'model_s');
end

[~, indModel, indNum] = intersect(cellfun(@str2num, model_u.genes), num(:, 1));
expressionData_u.gene(1:length(indModel)) = model_u.genes(indModel);
expressionData_u.value(1:length(indNum)) = num(indNum, 2);

[~, indModel, indNum] = intersect(cellfun(@str2num, model_c.genes), num(:, 1));
expressionData_c.gene(1:length(indModel)) = model_c.genes(indModel);
expressionData_c.value(1:length(indNum)) = num(indNum, 2);

[~, indModel, indNum] = intersect(cellfun(@str2num, model_s.genes), num(:, 1));
expressionData_s.gene(1:length(indModel)) = model_s.genes(indModel);
expressionData_s.value(1:length(indNum)) = num(indNum, 2);

if strcmp(figName,'U')
    %UNCONSTRAINED
    epsil = 1;
    tol = 1e-6;
    runtime = 7200;
    expressionCol = mapExpressionToReactions(model_u, expressionData_u);
    model = model_u;
end

if strcmp(figName,'C')
    % CONSTRAINED
    epsil = 1e-6;
    tol = 1e-8;
    runtime = 14400;
    expressionCol = mapExpressionToReactions(model_c, expressionData_c);
    model = model_c;    
end

if strcmp(figName,'S')
    %SEMI-CONSTRAINED
    epsil = 1;
    tol = 1e-6;
    runtime = 7200;
    expressionCol = mapExpressionToReactions(model_s, expressionData_s);
    model = model_s;
end

if strcmp(bb,'B')
    biomassRxn = model.rxns(strcmpi(model.rxns, 'biomass_reaction'));
    model = changeRxnBounds(model, biomassRxn, blb, 'l'); %Force biomass and ATP demand to be active
    atpDMInd = strncmp(model.rxns, 'DM_atp', 6) | strcmp(model.rxns, 'ATPM');
    figName = [figName,'B'];
end
if strcmp(bb,'F')
    model = addBiomassSinks(model);
    bmsInd=~cellfun(@isempty,strfind(model.rxns,'BMS_'));
    figName = [figName,'F'];
end
if strcmp(bb,'H')
    biomassRxn = model.rxns(strcmpi(model.rxns, 'biomass_reaction'));
    atpDMInd = strncmp(model.rxns, 'DM_atp', 6) | strcmp(model.rxns, 'ATPM');
    model = changeRxnBounds(model, biomassRxn, 1e-3, 'l'); %Force biomass and ATP demand to be active
    figName = [figName,'H'];
end

biomassRxnsInd = strncmpi(model.rxns, 'biomass', 7);
%w1
w1 = zeros(length(expressionCol),1);
w1(expressionCol >= 0) = 5*log(expressionCol(expressionCol>=0)/ths.p10);
w1(expressionCol < 0) = -2; % "unknown" entries get a weight of -2
w1(w1 < -max(w1)) = -max(w1);
if strcmp(bb,'B') || strcmp(bb,'H')
    w1(biomassRxnsInd | atpDMInd) = max(w1); %Biomass and ATP demand get high weight
else
    w1(bmsInd) = 0;
end

%w2
w2 = zeros(length(expressionCol),1);
w2(expressionCol >= 0) = 5*log(expressionCol(expressionCol>=0)/ths.mean);
w2(expressionCol < 0) = -2; % "unknown" entries get a weight of -2
w2(w2 < -max(w2)) = -max(w2);
if strcmp(bb,'B') || strcmp(bb,'H')
    w2(biomassRxnsInd | atpDMInd) = max(w2); %Biomass and ATP demand get high weight
else
    w2(bmsInd) = 0;
end

%w3
w3 = zeros(length(expressionCol),1);
w3(expressionCol >= 0) = 5*log(expressionCol(expressionCol>=0)/ths.p25);
w3(expressionCol < 0) = -2; % "unknown" entries get a weight of -2
w3(w3 < -max(w3)) = -max(w3);
if strcmp(bb,'B') || strcmp(bb,'H')
    w3(biomassRxnsInd | atpDMInd) = max(w3); %Biomass and ATP demand get high weight
else
    w3(bmsInd) = 0;
end

%w4
w4 = zeros(length(expressionCol),1);
w4(expressionCol >= 0) = 5*log(expressionCol(expressionCol>=0)/ths.p50);
w4(expressionCol < 0) = -2; % "unknown" entries get a weight of -2
w4(w4 < -max(w4)) = -max(w4);
if strcmp(bb,'B') || strcmp(bb,'H')
    w4(biomassRxnsInd | atpDMInd) = max(w4); %Biomass and ATP demand get high weight
else
    w4(bmsInd) = 0;
end

run_INIT(model, ths.p25, figName, epsil, w1, 1, modelName, tol, runtime, cellLine, overWrite);
run_INIT(model, ths.p50, figName, epsil, w2, 2, modelName, tol, runtime, cellLine, overWrite);
run_INIT(model, ths.p10, figName, epsil, w3, 3, modelName, tol, runtime, cellLine, overWrite);
run_INIT(model, ths.mean, figName, epsil, w4, 4, modelName, tol, runtime, cellLine, overWrite);
end

function run_INIT(model, ~, figName, epsil, w, id, modelName, tol, runtime, cellLine, overWrite)
    tName = ['INIT_', cellLine, '_', figName, num2str(id), '_', modelName];
    disp(tName)
    optionsLocal = struct('weights', {w}, 'tol', tol, 'runtime', runtime', ...
        'epsilon', epsil, 'solver', 'INIT', 'logfile', [tName,'.txt']);
    paramConsistency.epsilon=tol;
    paramConsistency.modeFlag=0;
    paramConsistency.method='fastcc';
    if (~overWrite && exist([pwd '/' tName '.mat'], 'file') ~=2)
        overWrite = 1;
    end
    if (overWrite)
        try
            cMod = createTissueSpecificModel(model, optionsLocal, 1, [], paramConsistency);
            cMod.name = tName;
            writeCbModel(cMod, 'mat', tName);
        catch ME
            warning('Failed to run INIT on model %s, figure %s with cell line %s', modelName, [figName num2str(id)], cellLine);
            warning(ME.message)
        end
    end
end