function runThresholdsmCADRE(figName, bb, modelName, cellLine)
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
        core = {};
        if strcmp(bb,'B')
            biomassRxnInd = strcmpi(model_u.rxns, 'biomass_reaction');
            biomassRxn = model_u.rxns(biomassRxnInd);
            atpDM = model_u.rxns(strncmp(model_u.rxns, 'DM_atp', 6) | strcmp(model_u.rxns, 'ATPM'));
            model_u = changeRxnBounds(model_u, biomassRxn, blb, 'l'); %Force biomass and ATP demand to be active
            core = [biomassRxn, atpDM];
            figName = [figName,'B'];
        end
        if strcmp(bb,'F')
            model_u = addBiomassSinks(model_u);
            figName = [figName,'F'];
            bmi = strcmpi(model_u.rxns, 'biomass_reaction');
            model_u.rxns{bmi} = 'Biomass_reaction_rename';
        end
        disp('UNCONSTRAINED MODEL')
        expressionCol_u = mapExpressionToReactions(model_u, expressionData_u);
        singleRun(core, model_u,  ths.p10, 1, 1/3, figName, 1, modelName, tol, expressionCol_u, cellLine)
        singleRun(core, model_u,  ths.mean, 1, 1/3, figName, 2, modelName, tol, expressionCol_u, cellLine)
        singleRun(core, model_u,  ths.p25, 1, 1/3, figName, 3, modelName, tol, expressionCol_u, cellLine)
        singleRun(core, model_u,  ths.p50, 1, 1/3, figName, 4, modelName, tol, expressionCol_u, cellLine)
    end
    if strcmp(figName,'C')
        tol = 1e-8;
        disp('CONSTRAINED MODEL')
        core = {};
        if strcmp(bb,'B')
            biomassRxnInd = strcmpi(model_c.rxns, 'biomass_reaction');
            biomassRxn = model_c.rxns(biomassRxnInd);
            atpDM = model_c.rxns(strncmp(model_c.rxns, 'DM_atp', 6) | strcmp(model_c.rxns, 'ATPM'));
            model_c = changeRxnBounds(model_c, biomassRxn, blb, 'l'); %Force biomass and ATP demand to be active
            core = [biomassRxn, atpDM];
            figName = [figName,'B'];
        end
        if strcmp(bb,'F')
            model_c = addBiomassSinks(model_c);
            figName = [figName,'F'];
            bmi = strcmpi(model_c.rxns, 'biomass_reaction');
            model_c.rxns{bmi} = 'Biomass_reaction_rename';
        end
        if strcmp(bb,'H')
            biomassRxnInd = strcmpi(model_c.rxns, 'biomass_reaction');
            biomassRxn = model_c.rxns{biomassRxnInd};
            atpDM = model_c.rxns(strncmp(model_c.rxns, 'DM_atp', 6) | strcmp(model_c.rxns, 'ATPM'));
            model_c = changeRxnBounds(model_c, biomassRxn, 1e-3, 'l'); %Force biomass and ATP demand to be active
            core = [biomassRxn, atpDM];            
            figName = [figName,'H'];
        end
        expressionCol_c = mapExpressionToReactions(model_c, expressionData_c);
        singleRun(core, model_c, ths.p10, 1, 1/3, figName, 1, modelName, tol, expressionCol_c, cellLine)
        singleRun(core, model_c, ths.mean, 1, 1/3, figName, 2, modelName, tol, expressionCol_c, cellLine)
        singleRun(core, model_c, ths.p25, 1, 1/3, figName, 3, modelName, tol, expressionCol_c, cellLine)
        singleRun(core, model_c, ths.p50, 1, 1/2, figName, 4, modelName, tol, expressionCol_c, cellLine)
    end
    if strcmp(figName,'S')
        tol = 1e-6;
        core = {};
        if strcmp(bb,'B')
            biomassRxnInd = strcmpi(model_s.rxns, 'biomass_reaction');
            biomassRxn = model_s.rxns(biomassRxnInd);
            atpDM = model_s.rxns(strncmp(model_s.rxns, 'DM_atp', 6) | strcmp(model_s.rxns, 'ATPM'));
            model_s = changeRxnBounds(model_s, biomassRxn, blb, 'l'); %Force biomass and ATP demand to be active
            core = [biomassRxn, atpDM];
            figName = [figName,'B'];
        end
        if strcmp(bb,'F')
            model_s = addBiomassSinks(model_s);
            figName = [figName,'F'];
            bmi = strcmpi(model_s.rxns, 'biomass_reaction');
            model_s.rxns{bmi} = 'Biomass_reaction_rename';
        end
        disp('SEMI-CONSTRAINED MODEL')
        expressionCol_s = mapExpressionToReactions(model_s, expressionData_s);
        singleRun(core, model_s, ths.p10, 1, 1/3, figName, 1, modelName, tol, expressionCol_s, cellLine)
        singleRun(core, model_s, ths.mean, 1, 1/3, figName, 2, modelName, tol, expressionCol_s, cellLine)
        singleRun(core, model_s, ths.p25, 1, 1/3, figName, 3, modelName, tol, expressionCol_s, cellLine)
        singleRun(core, model_s, ths.p50, 1, 1/3, figName, 4, modelName, tol, expressionCol_s, cellLine)
    end
    exit;
end

function singleRun(core, model, ht, mcheck, eta, figName, id, modelName, tol, expressionCol, cellLine)
    tName = ['mCADRE_',modelName, '_', cellLine, '_', figName, num2str(id)];
    disp(tName)
    disp('RUNNING mCADRE...')
     
    if (any(ismember(fieldnames(model), 'confidenceScores')))
        confidenceScores = str2double(model.confidenceScores);
    else
        confidenceScores = str2double(model.rxnConfidenceScores);
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
    try
        cMod = createTissueSpecificModel(model, optionsLocal); % consistency check not required, because mCADRE seems to run it on its own
        cMod.name = tName;
        writeCbModel(cMod, 'mat', tName);
    catch ME
        warning('Failed to run mCADRE on model %s, figure %s with cell line %s', modelName, [figName num2str(id)], cellLine);
        warning(ME.message)
    end
end