function runThresholdsMBA(figName, bb, modelName, cellLine)
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
    core = {};
    if strcmp(bb, 'B')
        biomassRxnInd = strcmpi(model_u.rxns, 'biomass_reaction');
        biomassRxn = model_u.rxns(biomassRxnInd);
        atpDM = model_u.rxns(strncmp(model_u.rxns, 'DM_atp', 6) | strcmp(model_u.rxns, 'ATPM'));
        model_u = changeRxnBounds(model_u, biomassRxn, blb, 'l'); %Force biomass and ATP demand to be active
        core = [biomassRxn, atpDM];
        figName = [figName,'B'];
    end
    if strcmp(bb, 'F')
        model_u = addBiomassSinks(model_u);
        figName = [figName,'F'];
    end
    expressionCol_u = mapExpressionToReactions(model_u, expressionData_u);
    run_MBA(core, model_u, expressionCol_u, figName, ths.p10, ths.p10, 1, modelName, tol, cellLine)
    run_MBA(core, model_u, expressionCol_u, figName, ths.mean, ths.mean, 2, modelName, tol, cellLine)
    run_MBA(core, model_u, expressionCol_u, figName, ths.p25, ths.p25, 3, modelName, tol, cellLine)
    run_MBA(core, model_u, expressionCol_u, figName, ths.p50, ths.p50, 4, modelName, tol, cellLine)
    run_MBA(core, model_u, expressionCol_u, figName, ths.mean, ths.p10, 5, modelName, tol, cellLine)
    run_MBA(core, model_u, expressionCol_u, figName, ths.p25, ths.p10, 6, modelName, tol, cellLine)
    run_MBA(core, model_u, expressionCol_u, figName, ths.p50, ths.p10, 7, modelName, tol, cellLine)
    run_MBA(core, model_u, expressionCol_u, figName, ths.p25, ths.mean, 8, modelName, tol, cellLine)
    run_MBA(core, model_u, expressionCol_u, figName, ths.p50, ths.mean, 9, modelName, tol, cellLine)
    run_MBA(core, model_u, expressionCol_u, figName, ths.p50, ths.p25, 10, modelName, tol, cellLine)
end

if strcmp(figName,'C')
    %CONSTRAINED
    tol = 1e-8;
    core = {};
    if strcmp(bb, 'B')
        biomassRxnInd = strcmpi(model_c.rxns, 'biomass_reaction');
        biomassRxn = model_c.rxns(biomassRxnInd);
        atpDM = model_c.rxns(strncmp(model_c.rxns, 'DM_atp', 6) | strcmp(model_c.rxns, 'ATPM'));
        model_c = changeRxnBounds(model_c, biomassRxn, blb, 'l'); %Force biomass and ATP demand to be active
        core = [biomassRxn, atpDM];
        figName = [figName,'B'];
    end
    if strcmp(bb, 'F')
        model_c = addBiomassSinks(model_c);
        figName = [figName,'F'];
    end
    if strcmp(bb, 'H')
        biomassRxnInd = strcmpi(model_c.rxns, 'biomass_reaction');
        biomassRxn = model_c.rxns(biomassRxnInd);
        atpDM = model_c.rxns(strncmp(model_c.rxns, 'DM_atp', 6) | strcmp(model_c.rxns, 'ATPM'));
        model_c = changeRxnBounds(model_c, biomassRxn, 1e-3, 'l'); %Force biomass and ATP demand to be active
        core = [biomassRxn, atpDM];
        figName = [figName,'H'];
    end
    expressionCol_c = mapExpressionToReactions(model_c, expressionData_c);
    run_MBA(core, model_c, expressionCol_c, figName, ths.p10, ths.p10, 1, modelName, tol, cellLine)
    run_MBA(core, model_c, expressionCol_c, figName, ths.mean, ths.mean, 2, modelName, tol, cellLine)
    run_MBA(core, model_c, expressionCol_c, figName, ths.p25, ths.p25, 3, modelName, tol, cellLine)
    run_MBA(core, model_c, expressionCol_c, figName, ths.p50, ths.p50, 4, modelName, tol, cellLine)
    run_MBA(core, model_c, expressionCol_c, figName, ths.mean, ths.p10, 5, modelName, tol, cellLine)
    run_MBA(core, model_c, expressionCol_c, figName, ths.p25, ths.p10, 6, modelName, tol, cellLine)
    run_MBA(core, model_c, expressionCol_c, figName, ths.p50, ths.p10, 7, modelName, tol, cellLine)
    run_MBA(core, model_c, expressionCol_c, figName, ths.p25, ths.mean, 8, modelName, tol, cellLine)
    run_MBA(core, model_c, expressionCol_c, figName, ths.p50, ths.mean, 9, modelName, tol, cellLine)
    run_MBA(core, model_c, expressionCol_c, figName, ths.p50, ths.p25, 10, modelName, tol, cellLine)
end
if strcmp(figName,'S')
    %SEMI-CONSTRAINED
    tol = 1e-6;
    core = {};
    if strcmp(bb, 'B')
        biomassRxnInd = strcmpi(model_s.rxns, 'biomass_reaction');
        biomassRxn = model_s.rxns(biomassRxnInd);
        atpDM = model_s.rxns(strncmp(model_s.rxns, 'DM_atp', 6) | strcmp(model_s.rxns, 'ATPM'));
        model_s = changeRxnBounds(model_s, biomassRxn, blb, 'l'); %Force biomass and ATP demand to be active
        core = [biomassRxn, atpDM];
        figName = [figName,'B'];
    end
    if strcmp(bb, 'F')
        model_s = addBiomassSinks(model_s);
        figName = [figName,'F'];
    end
    expressionCol_s = mapExpressionToReactions(model_s, expressionData_s);
    run_MBA(core, model_s, expressionCol_s, figName, ths.p10, ths.p10, 1, modelName, tol, cellLine)
    run_MBA(core, model_s, expressionCol_s, figName, ths.mean, ths.mean, 2, modelName, tol, cellLine)
    run_MBA(core, model_s, expressionCol_s, figName, ths.p25, ths.p25, 3, modelName, tol, cellLine)
    run_MBA(core, model_s, expressionCol_s, figName, ths.p50, ths.p50, 4, modelName, tol, cellLine)
    run_MBA(core, model_s, expressionCol_s, figName, ths.mean, ths.p10, 5, modelName, tol, cellLine)
    run_MBA(core, model_s, expressionCol_s, figName, ths.p25, ths.p10, 6, modelName, tol, cellLine)
    run_MBA(core, model_s, expressionCol_s, figName, ths.p50, ths.p10, 7, modelName, tol, cellLine)
    run_MBA(core, model_s, expressionCol_s, figName, ths.p25, ths.mean, 8, modelName, tol, cellLine)
    run_MBA(core, model_s, expressionCol_s, figName, ths.p50, ths.mean, 9, modelName, tol, cellLine)
    run_MBA(core, model_s, expressionCol_s, figName, ths.p50, ths.p25, 10, modelName, tol, cellLine)
end
exit;
end

function run_MBA(core, model, expressionCol, figName, mt, ut, id, modName, tol, cellLine)
    tName = ['MBA','_',figName, num2str(id),'_',modName];
    disp(tName)
    
    indH = find(expressionCol > ut);
    indM = find(expressionCol >= mt & expressionCol <= ut);
    CH = union(model.rxns(indH),core); %#ok<FNDSB>
    CM = model.rxns(indM); %#ok<FNDSB>
    NC = setdiff(model.rxns, union(CH, CM));
    %Biomass metabolite sinks are not pruned in call_MBA
    biomassMetSinksInd = contains(NC, 'BMS_'); 
    CM = union(CM, NC(biomassMetSinksInd));
    
    try
        cMod = MBA(model, CM, CH, tol);
        disp(['Number of rxns: ',num2str(numel(cMod.rxns))])
        cMod.name = tName;
        writeCbModel(cMod, 'mat', tName);
    catch ME
        warning('Failed to run MBA on model %s, figure %s with cell line %s', modelName, [figName num2str(id)], cellLine);
        warning(ME.message)
    end
        
end