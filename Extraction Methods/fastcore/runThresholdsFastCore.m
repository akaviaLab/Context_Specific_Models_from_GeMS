function runThresholdsFastCore(figName, bb, modName)
    %initCobraToolbox
        load(['ID_FPKM_', modName, '.mat'], 'num');
        load(['model_u_',modName,'.mat'], 'model_u')
    % This loads up the gene expression based on the FPKM file, so it has
    % the option of dealing with more reactions. There are three genes
    % (2878, 4259, 5337)
    % and four reactions (GTHPe, LTC4Sr, PCHOLPm_hs, PCHOLPr_hs) that are 
    % different based on the FPKM file to the existing gene_expr/gene_id structure.
    [~, indModel, indNum] = intersect(cellfun(@str2num, model_u.genes), num(:, 1));
    expressionData_u.gene(1:length(indModel)) = model_u.genes(indModel);
    expressionData_u.value(1:length(indNum)) = num(indNum, 2);
    % This sets up the original data
%      load(['gene_expr_u_',modName,'.mat'])
%      load(['gene_id_u_',modName,'.mat'])
%      expressionData_u2.gene = gene_id_u;
%      expressionData_u2.value = gene_expr_u;
    
    load(['model_c_',modName,'.mat'], 'model_c')
    [~, indModel, indNum] = intersect(cellfun(@str2num, model_c.genes), num(:, 1));
    expressionData_c.gene(1:length(indModel)) = model_c.genes(indModel);
    expressionData_c.value(1:length(indNum)) = num(indNum, 2);
    
    load(['model_s_',modName,'.mat'], 'model_s')
    [~, indModel, indNum] = intersect(cellfun(@str2num, model_s.genes), num(:, 1));
    expressionData_c.gene(1:length(indModel)) = model_s.genes(indModel);
    expressionData_c.value(1:length(indNum)) = num(indNum, 2);
    
    load(['growthRate_',modName,'.mat'], 'blb')
    load(['gene_threshold_',modName,'.mat'], 'ths')
    
    if strcmp(figName,'U')
        %UNCONSTRAINED
        core = {};
        if strcmp(bb,'B')
            model_u = changeRxnBounds(model_u, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
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
        singleRun(core, expressionCol, figName, model_u, epsil, ths.p10, 1, modName)
        singleRun(core, expressionCol, figName, model_u, epsil, ths.mean, 2, modName)
        singleRun(core, expressionCol, figName, model_u, epsil, ths.p25, 3, modName)
        singleRun(core, expressionCol, figName, model_u, epsil, ths.p50, 4, modName)
    end

    if strcmp(figName,'C')
        %CONSTRAINED
        core = {};
        if strcmp(bb,'B')
            model_c = changeRxnBounds(model_c, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
            core = {'Biomass_reaction','DM_atp(c)'};
            figName = [figName,'B'];
        end
        if strcmp(bb,'F')
            model_c = addBiomassSinks(model_c);
            figName = [figName,'F'];
        end
        if strcmp(bb,'H')
            model_c = changeRxnBounds(model_c, 'Biomass_reaction', 1e-3, 'l'); %Force biomass and ATP demand to be active
            core = {'Biomass_reaction','DM_atp(c)'};
            figName = [figName,'H'];
        end
        epsil = 1e-8;
        scaling = 1e3;
        expressionCol = mapExpressionToReactions(model_c, expressionData_c);
        singleRun(core, expressionCol, figName, model_c, epsil, ths.p10, 1, modName)
        singleRun(core, expressionCol, figName, model_c, epsil, ths.mean, 2, modName)
        singleRun(core, expressionCol, figName, model_c, epsil, ths.p25, 3, modName)
        singleRun(core, expressionCol, figName, model_c, epsil, ths.p50, 4, modName)
    end
    
    if strcmp(figName,'S')
        %CONSTRAINED
        core = {};
        if strcmp(bb,'B')
            model_s = changeRxnBounds(model_s, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
            core = {'Biomass_reaction','DM_atp(c)'};
            figName = [figName,'B'];
        end
        if strcmp(bb,'F')
            model_s = addBiomassSinks(model_s);
            figName = [figName,'F'];
        end
        epsil = 1e-6;
        expressionCol = mapExpressionToReactions(model_s, expressionData_s);
        singleRun(core, expressionCol, figName, model_s, epsil, ths.p10, 1, modName)
        singleRun(core, expressionCol, figName, model_s, epsil, ths.mean, 2, modName)
        singleRun(core, expressionCol, figName, model_s, epsil, ths.p25, 3, modName)
        singleRun(core, expressionCol, figName, model_s, epsil, ths.p50, 4, modName)
    end

end

function singleRun(core, expressionCol, figName, model, epsil, th, id, modName)
    tName = ['FastCore_',figName, num2str(id),'_',modName];
    disp(tName)
    C = find(expressionCol >= th);
    C = union(C, find(ismember(model.rxns, core)));
    cMod = fastcore(model, C, epsil);
    cMod.name = tName;
    writeCbModel(cMod, 'mat', tName);
end