function runThresholdsmCADRE(figName, bb, modelName, cellLine, overWrite)

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
    tol = 1e-6;
    disp('UNCONSTRAINED MODEL')
    expressionCol = mapExpressionToReactions(model_u, expressionData_u);
    model = model_u;
end
if strcmp(figName,'C')
    tol = 1e-8;
    disp('CONSTRAINED MODEL')
    expressionCol = mapExpressionToReactions(model_c, expressionData_c);
    model = model_c;
end
if strcmp(figName,'S')
    tol = 1e-6;
    disp('SEMI-CONSTRAINED MODEL')
    expressionCol = mapExpressionToReactions(model_s, expressionData_s);
    model = model_s;
end

core = {};
if strcmp(bb,'B')
    biomassRxnInd = strcmpi(model.rxns, 'biomass_reaction');
    biomassRxn = model.rxns(biomassRxnInd);
    atpDM = model.rxns(strncmp(model.rxns, 'DM_atp', 6) | strcmp(model.rxns, 'ATPM'));
    model = changeRxnBounds(model, biomassRxn, blb, 'l'); %Force biomass and ATP demand to be active
    core = [biomassRxn, atpDM];
    figName = [figName,'B'];
end
if strcmp(bb,'F')
    model = addBiomassSinks(model);
    figName = [figName,'F'];
    bmi = strcmpi(model.rxns, 'biomass_reaction');
    model.rxns{bmi} = 'Biomass_reaction_rename';
end
if strcmp(bb,'H')
    biomassRxnInd = strcmpi(model.rxns, 'biomass_reaction');
    biomassRxn = model.rxns{biomassRxnInd};
    atpDM = model.rxns(strncmp(model.rxns, 'DM_atp', 6) | strcmp(model.rxns, 'ATPM'));
    model = changeRxnBounds(model, biomassRxn, 1e-3, 'l'); %Force biomass and ATP demand to be active
    core = [biomassRxn, atpDM];
    figName = [figName,'H'];
end

singleRun(core, model, ths.p10, 1, 1/3, figName, 1, modelName, tol, expressionCol, cellLine, overWrite)
singleRun(core, model, ths.mean, 1, 1/3, figName, 2, modelName, tol, expressionCol, cellLine, overWrite)
singleRun(core, model, ths.p25, 1, 1/3, figName, 3, modelName, tol, expressionCol, cellLine, overWrite)
singleRun(core, model, ths.p50, 1, 1/2, figName, 4, modelName, tol, expressionCol, cellLine, overWrite)

end

function singleRun(core, model, ht, mcheck, eta, figName, id, modelName, tol, expressionCol, cellLine, overWrite)
    tName = ['mCADRE_', figName, num2str(id), '_', cellLine, '_', modelName];
    disp(tName)
    disp('RUNNING mCADRE...')
     
    if (any(ismember(fieldnames(model), 'confidenceScores')))
        confidenceScores = str2double(model.confidenceScores);
    else
        confidenceScores = model.rxnConfidenceScores;
    end
    confidenceScores(isnan(confidenceScores)) = 0;
    expressionCol = expressionCol / ht;
    expressionCol(expressionCol > 1) = 1;
    protectedRxns = union(core, model.rxns(expressionCol >= 1));
    optionsLocal = struct('solver', 'mCADRE',...
        'checkFunctionality', mcheck, 'eta', eta, 'tol', tol,...
        'ubiquityScore', {expressionCol},...
        'confidenceScores', {confidenceScores},...
        'protectedRxns', {protectedRxns});
    if (isempty(modelName))   % Recon1 can't produce the metabolites for the protected reaction
        optionsLocal.checkFunctionality = 0;
    end
    if (~overWrite && exist([pwd '/' tName '.mat'], 'file') ~=2)
        overWrite = 1;
    end
    if (overWrite)
        try
            cMod = createTissueSpecificModel(model, optionsLocal); % consistency check not required, because mCADRE seems to run it on its own
            cMod.name = tName;
            writeCbModel(cMod, 'mat', tName);
        catch ME
            warning('Failed to run mCADRE on model %s, figure %s with cell line %s', modelName, [figName num2str(id)], cellLine);
            warning(ME.message)
        end
    end
end