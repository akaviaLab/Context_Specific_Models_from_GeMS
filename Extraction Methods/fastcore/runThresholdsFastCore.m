function runThresholdsFastCore(figName, bb, modelName, cellLine, overWrite)

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
    epsil = 1e-6;
    scaling = 1e3;
    expressionCol = mapExpressionToReactions(model_u, expressionData_u);
    model = model_u;
end

if strcmp(figName,'C')
    %CONSTRAINED  
    epsil = 1e-8;
    scaling = 1e3;
    expressionCol = mapExpressionToReactions(model_c, expressionData_c);
    model = model_c;
end

if strcmp(figName,'S')
    %SEMI-CONSTRAINED
    epsil = 1e-6;
    scaling = 1e3;
    expressionCol = mapExpressionToReactions(model_s, expressionData_s);
    model = model_s;
end

core = [];
if strcmp(bb,'B')
    biomassRxn = model.rxns(strcmpi(model.rxns, 'biomass_reaction'));
    model = changeRxnBounds(model, biomassRxn, blb, 'l'); %Force biomass and ATP demand to be active
    atpDMInd = strncmp(model.rxns, 'DM_atp', 6) | strcmp(model.rxns, 'ATPM');
    biomassRxnsInd = strncmpi(model.rxns, 'biomass', 7);
    core = find(biomassRxnsInd | atpDMInd);
    figName = [figName,'B'];
end
if strcmp(bb,'F')
    model = addBiomassSinks(model);
    figName = [figName,'F'];
end
if strcmp(bb,'H')
    biomassRxn = model.rxns(strcmpi(model.rxns, 'biomass_reaction'));
    model = changeRxnBounds(model, biomassRxn, 1e-3, 'l'); %Force biomass and ATP demand to be active
    biomassRxnsInd = strncmpi(model.rxns, 'biomass', 7);
    atpDMInd = strncmp(model.rxns, 'DM_atp', 6) | strcmp(model.rxns, 'ATPM');
    core = find(biomassRxnsInd | atpDMInd);
    figName = [figName,'H'];
end

run_fastcore(core, expressionCol, figName, model, epsil, scaling, ths.p10, 1, modelName, cellLine, overWrite)
run_fastcore(core, expressionCol, figName, model, epsil, scaling, ths.mean, 2, modelName, cellLine, overWrite)
run_fastcore(core, expressionCol, figName, model, epsil, scaling, ths.p25, 3, modelName, cellLine, overWrite)
run_fastcore(core, expressionCol, figName, model, epsil, scaling, ths.p50, 4, modelName, cellLine, overWrite)

end

function run_fastcore(core, expressionCol, figName, model, epsil, scaling, th, id, modelName, cellLine, overWrite)
tName = ['fastcore_', cellLine, '_', figName, num2str(id), '_', modelName];
disp(tName)
C = find(expressionCol >= th);
C = union(C, core);
optionsLocal = struct('solver', 'fastCore', 'core', {C}, 'epsilon', epsil, 'printLevel', 2);
if (~overWrite && exist([pwd '/' tName '.mat'], 'file') ~=2)
    overWrite = 1;
end
if (overWrite)
    try
        cMod = createTissueSpecificModel(model, optionsLocal);
        cMod.name = tName;
        writeCbModel(cMod, 'mat', tName);
    catch ME
        warning('Failed fastcore on model %s, figure %s with cell line %s', modelName, [figName num2str(id)], cellLine);
        warning(ME.message)
    end
end
end