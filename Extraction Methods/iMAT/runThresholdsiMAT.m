function runThresholdsiMAT(figName, bb, modelName, cellLine, overWrite)

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
    tol = 1e-6;
    runtime = 3600;
    epsil = 1;
    expressionCol = mapExpressionToReactions(model_u, expressionData_u);
    model = model_u;
end

if strcmp(figName,'C')
    %CONSTRAINED
    tol = 1e-8;
    runtime = 14400; % Works with 32 CPUs, 96G of memory, otherwise, can take longer
    epsil = 1e-6;
    expressionCol = mapExpressionToReactions(model_c, expressionData_c);
    model = model_c;
end

if strcmp(figName,'S')
    %SEMI-CONSTRAINED
    tol = 1e-6;
    runtime = 7200;
    epsil = 1;
    expressionCol = mapExpressionToReactions(model_s, expressionData_s);
    model = model_s;
end

if strcmp(bb, 'B')
    biomassRxn = model.rxns(strcmpi(model.rxns, 'biomass_reaction'));
    model = changeRxnBounds(model, biomassRxn, blb, 'l'); %Force biomass and ATP demand to be active
    biomassRxnsInd = strncmpi(model.rxns, 'biomass', 7);
    atpDMInd = strncmp(model.rxns, 'DM_atp', 6) | strcmp(model.rxns, 'ATPM');
    core = model.rxns(biomassRxnsInd | atpDMInd);
    figName = [figName,'B'];
end
if strcmp(bb, 'F')
    model = addBiomassSinks(model);
    figName = [figName,'F'];
end
if strcmp(bb, 'H')
    biomassRxn = model.rxns(strcmpi(model.rxns, 'biomass_reaction'));
    model = changeRxnBounds(model, biomassRxn, 1e-3, 'l'); %Force biomass and ATP demand to be active
    biomassRxnsInd = strncmpi(model.rxns, 'biomass', 7);
    atpDMInd = strncmp(model.rxns, 'DM_atp', 6) | strcmp(model.rxns, 'ATPM');
    core = model.rxns(biomassRxnsInd | atpDMInd);
    figName = [figName,'H'];
end

run_iMat(core, model, expressionCol, figName, epsil, ths.p10, ths.p10, 1, modelName, tol, runtime, cellLine, overWrite)
run_iMat(core, model, expressionCol, figName, epsil, ths.mean, ths.mean, 2, modelName, tol, runtime, cellLine, overWrite)
run_iMat(core, model, expressionCol, figName, epsil, ths.p25, ths.p25, 3, modelName, tol, runtime, cellLine, overWrite)
run_iMat(core, model, expressionCol, figName, epsil, ths.p50, ths.p50, 4, modelName, tol, runtime, cellLine, overWrite)
run_iMat(core, model, expressionCol, figName, epsil, ths.mean, ths.p10, 5, modelName, tol, runtime, cellLine, overWrite)
run_iMat(core, model, expressionCol, figName, epsil, ths.p25, ths.p10, 6, modelName, tol, runtime, cellLine, overWrite)
run_iMat(core, model, expressionCol, figName, epsil, ths.p50, ths.p10, 7, modelName, tol, runtime, cellLine, overWrite)
run_iMat(core, model, expressionCol, figName, epsil, ths.p25, ths.mean, 8, modelName, tol, runtime, cellLine, overWrite)
run_iMat(core, model, expressionCol, figName, epsil, ths.p50, ths.mean, 9, modelName, tol, runtime, cellLine, overWrite)
run_iMat(core, model, expressionCol, figName, epsil, ths.p50, ths.p25, 10, modelName, tol, runtime, cellLine, overWrite)
end

function run_iMat(core, model, expressionCol, figName, epsil, lb, ub, id, modelName, tol, runtime, cellLine, overWrite)
    tName = ['iMAT_', cellLine, '_', figName, num2str(id), '_', modelName];
    disp(tName)
    paramConsistency.epsilon=tol;
    paramConsistency.modeFlag=0;
    paramConsistency.method='fastcc';
    optionsLocal = struct('solver', 'iMAT', 'expressionRxns', {expressionCol}, 'threshold_lb', lb, 'epsilon', epsil,...
        'threshold_ub', ub, 'tol', tol, 'core', {core}, 'logfile', [tName,'_2.txt'], 'runtime', runtime);
    if (~overWrite && exist([pwd '/' tName '.mat'], 'file') ~=2)
        overWrite = 1;
    end
    if (overWrite)
        try
            tic
            cMod2 = createTissueSpecificModel(model, optionsLocal, 1, [], paramConsistency);
            toc
            cMod2.name = tName;
            writeCbModel(cMod2, 'mat', tName)
        catch ME
            warning('Failed to run iMAT on model %s, figure %s with cell line %s', modelName, [figName num2str(id)], cellLine);
            warning(ME.message)
        end
    end
end