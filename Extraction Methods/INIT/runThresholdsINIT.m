function runThresholdsINIT(figName, bb, modelName, cellLine)
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
    tol = 1e-6;
    runtime = 3600;
    if strcmp(bb,'B')
        biomassRxnInd = strcmpi(model_u.rxns, 'biomass_reaction');
        biomassRxn = model_u.rxns(biomassRxnInd);
        atpDM = model_u.rxns(strncmp(model_u.rxns, 'DM_atp', 6) | strcmp(model_u.rxns, 'ATPM'));
        model = changeRxnBounds(model_u, biomassRxn, blb, 'l'); %Force biomass and ATP demand to be active
        figName = [figName,'B'];
    end
    if strcmp(bb,'F')
        model = addBiomassSinks(model_u);
        bmsInd=~cellfun(@isempty,strfind(model.rxns,'BMS_'));
        figName = [figName,'F'];
    end
    expressionCol = mapExpressionToReactions(model_u, expressionData_u);
end

if strcmp(figName,'C')
    tol = 1e-8;
    runtime = 7200;
    if strcmp(bb,'B')
        biomassRxnInd = strcmpi(model_c.rxns, 'biomass_reaction');
        biomassRxn = model_c.rxns(biomassRxnInd);
        atpDM = model_c.rxns(strncmp(model_c.rxns, 'DM_atp', 6) | strcmp(model_c.rxns, 'ATPM'));
        model = changeRxnBounds(model_c, biomassRxn, blb, 'l'); %Force biomass and ATP demand to be active
        figName = [figName,'B'];
    end
    if strcmp(bb,'F')
        model = addBiomassSinks(model_c);
        bmsInd=~cellfun(@isempty,strfind(model.rxns,'BMS_'));
        figName = [figName,'F'];
    end
    if strcmp(bb,'H')
        biomassRxnInd = strcmpi(model_c.rxns, 'biomass_reaction');
        biomassRxn = model_c.rxns(biomassRxnInd);
        atpDM = model_c.rxns(strncmp(model_c.rxns, 'DM_atp', 6) | strcmp(model_c.rxns, 'ATPM'));
        model = changeRxnBounds(model_c, biomassRxn, 1e-3, 'l'); %Force biomass and ATP demand to be active
        figName = [figName,'H'];
    end
    expressionCol = mapExpressionToReactions(model_c, expressionData_c);
end

if strcmp(figName,'S')
    tol = 1e-6;
    runtime = 3600;
    if strcmp(bb,'B')
        biomassRxnInd = strcmpi(model_s.rxns, 'biomass_reaction');
        biomassRxn = model_s.rxns(biomassRxnInd);
        atpDM = model_s.rxns(strncmp(model_s.rxns, 'DM_atp', 6) | strcmp(model_s.rxns, 'ATPM'));
        model = changeRxnBounds(model_s, biomassRxn, blb, 'l'); %Force biomass and ATP demand to be active
        figName = [figName,'B'];
    end
    if strcmp(bb,'F')
        model = addBiomassSinks(model_s);
        bmsInd=~cellfun(@isempty,strfind(model.rxns,'BMS_'));
        figName = [figName,'F'];
    end
    expressionCol = mapExpressionToReactions(model_s, expressionData_s);
end


%w1
w1 = zeros(length(expressionCol),1);
w1(expressionCol >= 0) = 5*log(expressionCol(expressionCol>=0)/ths.p10);
w1(expressionCol < 0) = -2; % "unknown" entries get a weight of -2
w1(w1 < -max(w1)) = -max(w1);
if strcmp(bb,'B') || strcmp(bb,'H')
    w1(ismember(model.rxns, [biomassRxn, atpDM])) = max(w1); %Biomass and ATP demand get high weight
else
    w1(bmsInd) = 0;
end

%w2
w2 = zeros(length(expressionCol),1);
w2(expressionCol >= 0) = 5*log(expressionCol(expressionCol>=0)/ths.mean);
w2(expressionCol < 0) = -2; % "unknown" entries get a weight of -2
w2(w2 < -max(w2)) = -max(w2);
if strcmp(bb,'B') || strcmp(bb,'H')
    w2(ismember(model.rxns, [biomassRxn, atpDM])) = max(w2); %Biomass and ATP demand get high weight
end

%w3
w3 = zeros(length(expressionCol),1);
w3(expressionCol >= 0) = 5*log(expressionCol(expressionCol>=0)/ths.p25);
w3(expressionCol < 0) = -2; % "unknown" entries get a weight of -2
w3(w3 < -max(w3)) = -max(w3);
if strcmp(bb,'B') || strcmp(bb,'H')
    w3(ismember(model.rxns, [biomassRxn, atpDM])) = max(w3); %Biomass and ATP demand get high weight
else
    w3(bmsInd) = 0;
end

%w4
w4 = zeros(length(expressionCol),1);
w4(expressionCol >= 0) = 5*log(expressionCol(expressionCol>=0)/ths.p50);
w4(expressionCol < 0) = -2; % "unknown" entries get a weight of -2
w4(w4 < -max(w4)) = -max(w4);
if strcmp(bb,'B') || strcmp(bb,'H')
    w4(ismember(model.rxns, [biomassRxn, atpDM])) = max(w4); %Biomass and ATP demand get high weight
else
    w4(bmsInd) = 0;
end
    
if strcmp(figName,'UB') || strcmp(figName,'UF')
    %UNCONSTRAINED
    epsil = 1;
    run_INIT(model,ths.p25, figName, epsil, w1, 1, modelName, tol, runtime, cellLine);
    run_INIT(model,ths.p50, figName, epsil, w2, 2, modelName, tol, runtime, cellLine);
    run_INIT(model,ths.p10, figName, epsil, w3, 3, modelName, tol, runtime, cellLine);
    run_INIT(model,ths.mean, figName, epsil, w4, 4, modelName, tol, runtime, cellLine);
end

if strcmp(figName,'CB') || strcmp(figName,'CF') || strcmp(figName,'CH')
    %CONSTRAINED
    epsil = 1e-6;
    run_INIT(model,ths.p25, figName, epsil, w1, 1, modelName, tol, runtime, cellLine);
    run_INIT(model,ths.p50, figName, epsil, w2, 2, modelName, tol, runtime, cellLine);
    run_INIT(model,ths.p10, figName, epsil, w3, 3, modelName, tol, runtime, cellLine);
    run_INIT(model,ths.mean, figName, epsil, w4, 4, modelName, tol, runtime, cellLine);
end
if strcmp(figName,'SB') || strcmp(figName,'SF')
    %SEMI-CONSTRAINED
    epsil = 1;
    run_INIT(model,ths.p25, figName, epsil, w1, 1, modelName, tol, runtime, cellLine);
    run_INIT(model,ths.p50, figName, epsil, w2, 2, modelName, tol, runtime, cellLine);
    run_INIT(model,ths.p10, figName, epsil, w3, 3, modelName, tol, runtime, cellLine);
    run_INIT(model,ths.mean, figName, epsil, w4, 4, modelName, tol, runtime, cellLine);

end
exit;
end

function run_INIT(model, ~, figName, epsil, w, id, modelName, tol, runtime, cellLine)
    tName = ['INIT_',figName, num2str(id),'_',modelName,'_',cellLine];
    disp(tName)
    cMod = call_INIT(model, epsil, w, tol, [tName,'.txt'], runtime);
    cMod.name = tName;
    writeCbModel(cMod, 'mat', tName);
    optionsLocal = struct('weights', {w}, 'tol', tol, 'runtime', runtime', ...
        'solver', 'INIT', 'logfile', [tName,'.txt']);
    paramConsistency.epsilon=tol;
    paramConsistency.modeFlag=0;
    paramConsistency.method='fastcc';
    try
        cMod2 = createTissueSpecificModel(model, optionsLocal, 1, [], paramConsistency);
        cMod2.name = tName;
        writeCbModel(cMod2, 'mat', [tName '_2']);
        if (~isSameCobraModel(cMod, cMod2))
            fprintf('When running INIT with model %s, fig %s and cell line %s, the old and new models are different!\n', modelName, figName, cellLine);
        end
    catch ME
        warning('Failed to run INIT on model %s, figure %s with cell line %s', modelName, [figName num2str(id)], cellLine);
        warning(ME.message)
    end
end