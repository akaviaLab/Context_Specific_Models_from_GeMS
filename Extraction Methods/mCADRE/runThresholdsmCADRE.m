function runThresholdsmCADRE(figName, bb, modName)
initCobraToolbox
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
            model_u = changeRxnBounds(model_u, biomassRxn, blb, 'l'); %Force biomass and ATP demand to be active
            core = {biomassRxn,'DM_atp(c)'};
            figName = [figName,'B'];
        end
        if strcmp(bb,'F')
            model_u = addBiomassSinks(model_u);
            figName = [figName,'F'];
            bmi = ismember(model_u.rxns,'Biomass_reaction');
            model_u.rxns{bmi} = 'Biomass_reaction_rename';
        end
        disp('UNCONSTRAINED MODEL')
        expressionCol_u = mapExpressionToReactions(model_u, expressionData_u);
        singleRun(core, model_u, gene_id_u, gene_expr_u, parsedGPR_u, corrRxn_u, ths.p10, 1, 1/3, figName, 1, modName, tol)
        singleRun(core, model_u, gene_id_u, gene_expr_u, parsedGPR_u, corrRxn_u, ths.mean, 1, 1/3, figName, 2, modName, tol)
        singleRun(core, model_u, gene_id_u, gene_expr_u, parsedGPR_u, corrRxn_u, ths.p25, 1, 1/3, figName, 3, modName, tol)
        singleRun(core, model_u, gene_id_u, gene_expr_u, parsedGPR_u, corrRxn_u, ths.p50, 1, 1/3, figName, 4, modName, tol)
    end
    if strcmp(figName,'C')
        tol = 1e-8;
        disp('CONSTRAINED MODEL')
        core = {};
        if strcmp(bb,'B')
            biomassRxnInd = strcmpi(model_c.rxns, 'biomass_reaction');
            biomassRxn = model_c.rxns(biomassRxnInd);
            model_c = changeRxnBounds(model_c, biomassRxn, blb, 'l'); %Force biomass and ATP demand to be active
            core = {biomassRxn,'DM_atp(c)'};
            figName = [figName,'B'];
        end
        if strcmp(bb,'F')
            model_c = addBiomassSinks(model_c);
            figName = [figName,'F'];
            bmi = ismember(model_c.rxns,'Biomass_reaction');
            model_c.rxns{bmi} = 'Biomass_reaction_rename';
        end
        if strcmp(bb,'H')
            biomassRxnInd = strcmpi(model_c.rxns, 'biomass_reaction');
            biomassRxn = model_c.rxns(biomassRxnInd);
            model_c = changeRxnBounds(model_c, biomassRxn, blb, 'l'); %Force biomass and ATP demand to be active
            core = {biomassRxn,'DM_atp(c)'};
            figName = [figName,'H'];
        end
        expressionCol_c = mapExpressionToReactions(model_c, expressionData_c);
        singleRun(core, model_c, gene_id_c, gene_expr_c, parsedGPR_c, corrRxn_c, ths.p10, 1, 1/3, figName, 1, modName, tol)
        singleRun(core, model_c, gene_id_c, gene_expr_c, parsedGPR_c, corrRxn_c, ths.mean, 1, 1/3, figName, 2, modName, tol)
        singleRun(core, model_c, gene_id_c, gene_expr_c, parsedGPR_c, corrRxn_c, ths.p25, 1, 1/3, figName, 3, modName, tol)
        singleRun(core, model_c, gene_id_c, gene_expr_c, parsedGPR_c, corrRxn_c, ths.p50, 1, 1/2, figName, 4, modName, tol)
    end
    if strcmp(figName,'S')
        tol = 1e-6;
        core = {};
        if strcmp(bb,'B')
            biomassRxnInd = strcmpi(model_s.rxns, 'biomass_reaction');
            biomassRxn = model_s.rxns(biomassRxnInd);
            model_s = changeRxnBounds(model_s, biomassRxn, blb, 'l'); %Force biomass and ATP demand to be active
            core = {biomassRxn,'DM_atp(c)'};
            figName = [figName,'B'];
        end
        if strcmp(bb,'F')
            model_s = addBiomassSinks(model_s);
            figName = [figName,'F'];
            bmi = ismember(model_s.rxns,'Biomass_reaction');
            model_s.rxns{bmi} = 'Biomass_reaction_rename';
        end
        disp('SEMI-CONSTRAINED MODEL')
        expressionCol_s = mapExpressionToReactions(model_s, expressionData_s);
        singleRun(core, model_s, gene_id_s, gene_expr_s, parsedGPR_s, corrRxn_s, ths.p10, 1, 1/3, figName, 1, modName, tol)
        singleRun(core, model_s, gene_id_s, gene_expr_s, parsedGPR_s, corrRxn_s, ths.mean, 1, 1/3, figName, 2, modName, tol)
        singleRun(core, model_s, gene_id_s, gene_expr_s, parsedGPR_s, corrRxn_s, ths.p25, 1, 1/3, figName, 3, modName, tol)
        singleRun(core, model_s, gene_id_s, gene_expr_s, parsedGPR_s, corrRxn_s, ths.p50, 1, 1/3, figName, 4, modName, tol)
    end
    exit;
end

function singleRun(core, model, gene_names, gene_exp, parsedGPR, corrRxn, ht, mcheck, eta, figName, id, modName, tol)
    tName = ['mCADRE_',figName,num2str(id),'_',modName];
    disp(tName)
    disp('RUNNING mCADRE...')
    cMod=call_mCADRE(model, gene_names, gene_exp, parsedGPR, corrRxn, core, ht, mcheck, eta, tol);
    inactiveRxns = CheckModelConsistency(cMod, tol);
    %Make sure output model is consistent
    if ~isempty(inactiveRxns)
        errmsg = 'Output model is inconsistent';
        save(['INC_',tName,'.mat'],'errmsg')
    end
    cMod.name = tName;
    writeCbModel([tName,'.mat'],tName)
end