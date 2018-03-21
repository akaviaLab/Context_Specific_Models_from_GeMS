function runThresholdsFastCore(figName, bb, modelName, cellLine)
    %initCobraToolbox
    % This loads up the gene expression based on the FPKM file, so it has
    % the option of dealing with more reactions. There are three genes
    % (2878, 4259, 5337)
    % and four reactions (GTHPe, LTC4Sr, PCHOLPm_hs, PCHOLPr_hs) that are 
    % different based on the FPKM file to the existing gene_expr/gene_id structure.
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
    expressionData_c.gene(1:length(indModel)) = model_s.genes(indModel);
    expressionData_c.value(1:length(indNum)) = num(indNum, 2);
        
    if strcmp(figName,'U')
        %UNCONSTRAINED
        core = {};
        if strcmp(bb,'B')
            biomassRxn = strcmpi(model_u.rxns, 'biomass_reaction');
            model_u = changeRxnBounds(model_u, model_u.rxns(biomassRxn), blb, 'l'); %Force biomass and ATP demand to be active
            core = {'Biomass_reaction','DM_atp(c)'};
            figName = [figName,'B'];
        end
        if strcmp(bb,'F')
            model_u = addBiomassSinks(model_u);
            figName = [figName,'F'];
        end
        epsil = 1e-6;
        scaling = 1e3;
        expressionCol = mapExpressionToReactions(model_u, expressionData_u);
        singleRun(core, expressionCol, figName, model_u, epsil, ths.p10, 1, modelName, cellLine)
        singleRun(core, expressionCol, figName, model_u, epsil, ths.mean, 2, modelName, cellLine)
        singleRun(core, expressionCol, figName, model_u, epsil, ths.p25, 3, modelName, cellLine)
        singleRun(core, expressionCol, figName, model_u, epsil, ths.p50, 4, modelName, cellLine)
    end

    if strcmp(figName,'C')
        %CONSTRAINED
        core = {};
        biomassRxn = strcmpi(model_c.rxns, 'biomass_reaction');
        if strcmp(bb,'B')
            model_c = changeRxnBounds(model_c, model_c.rxns(biomassRxn), blb, 'l'); %Force biomass and ATP demand to be active
            core = {'Biomass_reaction','DM_atp(c)'};
            figName = [figName,'B'];
        end
        if strcmp(bb,'F')
            model_c = addBiomassSinks(model_c);
            figName = [figName,'F'];
        end
        if strcmp(bb,'H')
            model_c = changeRxnBounds(model_c, model_c.rxns(biomassRxn), 1e-3, 'l'); %Force biomass and ATP demand to be active
            core = {'Biomass_reaction','DM_atp(c)'};
            figName = [figName,'H'];
        end
        epsil = 1e-8;
        scaling = 1e3;
        expressionCol = mapExpressionToReactions(model_c, expressionData_c);
        singleRun(core, expressionCol, figName, model_c, epsil, ths.p10, 1, modelName, cellLine)
        singleRun(core, expressionCol, figName, model_c, epsil, ths.mean, 2, modelName, cellLine)
        singleRun(core, expressionCol, figName, model_c, epsil, ths.p25, 3, modelName, cellLine)
        singleRun(core, expressionCol, figName, model_c, epsil, ths.p50, 4, modelName, cellLine)
    end
    
    if strcmp(figName,'S')
        %CONSTRAINED
        core = {};
        if strcmp(bb,'B')
            biomassRxn = strcmpi(model_s.rxns, 'biomass_reaction');
            model_s = changeRxnBounds(model_s, model_s.rxns(biomassRxn), blb, 'l'); %Force biomass and ATP demand to be active
            core = {'Biomass_reaction','DM_atp(c)'};
            figName = [figName,'B'];
        end
        if strcmp(bb,'F')
            model_s = addBiomassSinks(model_s);
            figName = [figName,'F'];
        end
        epsil = 1e-6;
        expressionCol = mapExpressionToReactions(model_s, expressionData_s);
        singleRun(core, expressionCol, figName, model_s, epsil, ths.p10, 1, modelName, cellLine)
        singleRun(core, expressionCol, figName, model_s, epsil, ths.mean, 2, modelName, cellLine)
        singleRun(core, expressionCol, figName, model_s, epsil, ths.p25, 3, modelName, cellLine)
        singleRun(core, expressionCol, figName, model_s, epsil, ths.p50, 4, modelName, cellLine)
    end

end

function singleRun(core, expressionCol, figName, model, epsil, th, id, modelName, cellLine)
    tName = ['FastCore_',modelName, '_', figName, num2str(id),'_',cellLine];
    disp(tName)
    C = find(expressionCol >= th);
    C = union(C, find(ismember(model.rxns, core)));
    try
    cMod = fastcore(model, C, epsil);
    cMod.name = tName;
    writeCbModel(cMod, 'mat', tName);
    catch ME
        warning('Failed to run fastcore on model %s, figure %s with cell line %s', modelName, [figName num2str(id)], cellLine);
        warning(ME.message)
    end
end