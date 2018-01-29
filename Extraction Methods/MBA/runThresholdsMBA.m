function runThresholdsMBA(figName, bb, modName)
load(['ID_FPKM_', modName, '.mat'], 'num');
load(['model_u_',modName,'.mat'], 'model_u')
[~, indModel, indNum] = intersect(cellfun(@str2num, model_u.genes), num(:, 1));
expressionData_u.gene(1:length(indModel)) = model_u.genes(indModel);
expressionData_u.value(1:length(indNum)) = num(indNum, 2);

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
    tol = 1e-6;
    core = {};
    if strcmp(bb, 'B')
        model_u = changeRxnBounds(model_u, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
        core = {'Biomass_reaction','DM_atp(c)'};
        figName = [figName,'B'];
    end
    if strcmp(bb, 'F')
        model_u = addBiomassSinks(model_u);
        figName = [figName,'F'];
    end
     expressionCol = mapExpressionToReactions(model_u, expressionData_u);
    run_MBA(core, model_u, expressionCol, figName, ths.p10, ths.p10, 1, modName, tol)
    run_MBA(core, model_u, expressionCol, figName, ths.mean, ths.mean, 2, modName, tol)
    run_MBA(core, model_u, expressionCol, figName, ths.p25, ths.p25, 3, modName, tol)
    run_MBA(core, model_u, expressionCol, figName, ths.p50, ths.p50, 4, modName, tol)
    run_MBA(core, model_u, expressionCol, figName, ths.mean, ths.p10, 5, modName, tol)
    run_MBA(core, model_u, expressionCol, figName, ths.p25, ths.p10, 6, modName, tol)
    run_MBA(core, model_u, expressionCol, figName, ths.p50, ths.p10, 7, modName, tol)
    run_MBA(core, model_u, expressionCol, figName, ths.p25, ths.mean, 8, modName, tol)
    run_MBA(core, model_u, expressionCol, figName, ths.p50, ths.mean, 9, modName, tol)
    run_MBA(core, model_u, expressionCol, figName, ths.p50, ths.p25, 10, modName, tol)
end

if strcmp(figName,'C')
    %CONSTRAINED
    tol = 1e-8;
    core = {};
    if strcmp(bb, 'B')
        model_c = changeRxnBounds(model_c, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
        core = {'Biomass_reaction','DM_atp(c)'};
        figName = [figName,'B'];
    end
    if strcmp(bb, 'F')
        model_c = addBiomassSinks(model_c);
        figName = [figName,'F'];
    end
    if strcmp(bb, 'H')
        model_c = changeRxnBounds(model_c, 'Biomass_reaction', 1e-3, 'l'); %Force biomass and ATP demand to be active
        core = {'Biomass_reaction','DM_atp(c)'};
        figName = [figName,'H'];
    end
    expressionCol = mapExpressionToReactions(model_c, expressionData_c);
    run_MBA(core, model_c, expressionCol, figName, ths.p10, ths.p10, 1, modName, tol)
    run_MBA(core, model_c, expressionCol, figName, ths.mean, ths.mean, 2, modName, tol)
    run_MBA(core, model_c, expressionCol, figName, ths.p25, ths.p25, 3, modName, tol)
    run_MBA(core, model_c, expressionCol, figName, ths.p50, ths.p50, 4, modName, tol)
    run_MBA(core, model_c, expressionCol, figName, ths.mean, ths.p10, 5, modName, tol)
    run_MBA(core, model_c, expressionCol, figName, ths.p25, ths.p10, 6, modName, tol)
    run_MBA(core, model_c, expressionCol, figName, ths.p50, ths.p10, 7, modName, tol)
    run_MBA(core, model_c, expressionCol, figName, ths.p25, ths.mean, 8, modName, tol)
    run_MBA(core, model_c, expressionCol, figName, ths.p50, ths.mean, 9, modName, tol)
    run_MBA(core, model_c, expressionCol, figName, ths.p50, ths.p25, 10, modName, tol)
end
if strcmp(figName,'S')
    %SEMI-CONSTRAINED
    tol = 1e-6;
    core = {};
    if strcmp(bb, 'B')
        model_s = changeRxnBounds(model_s, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
        core = {'Biomass_reaction','DM_atp(c)'};
        figName = [figName,'B'];
    end
    if strcmp(bb, 'F')
        model_s = addBiomassSinks(model_s);
        figName = [figName,'F'];
    end
     expressionCol = mapExpressionToReactions(model_s, expressionData_s);
    run_MBA(core, model_s, expressionCol, figName, ths.p10, ths.p10, 1, modName, tol)
    run_MBA(core, model_s, expressionCol, figName, ths.mean, ths.mean, 2, modName, tol)
    run_MBA(core, model_s, expressionCol, figName, ths.p25, ths.p25, 3, modName, tol)
    run_MBA(core, model_s, expressionCol, figName, ths.p50, ths.p50, 4, modName, tol)
    run_MBA(core, model_s, expressionCol, figName, ths.mean, ths.p10, 5, modName, tol)
    run_MBA(core, model_s, expressionCol, figName, ths.p25, ths.p10, 6, modName, tol)
    run_MBA(core, model_s, expressionCol, figName, ths.p50, ths.p10, 7, modName, tol)
    run_MBA(core, model_s, expressionCol, figName, ths.p25, ths.mean, 8, modName, tol)
    run_MBA(core, model_s, expressionCol, figName, ths.p50, ths.mean, 9, modName, tol)
    run_MBA(core, model_s, expressionCol, figName, ths.p50, ths.p25, 10, modName, tol)
end
exit;
end

function run_MBA(core, model, expressionCol, figName, mt, ut, id, modName, tol)
    paramConsistency.epsilon=1e-10;
    paramConsistency.modeFlag=0;
    paramConsistency.method='fastcc';
    tName = ['MBA','_',figName, num2str(id),'_',modName];
    disp(tName)
    %Make sure output model is consistent, if not, run again
    incon = true;
    nrun = 0;
    while incon && nrun < 10
            indH = find(expressionCol > ut);
    indM = find(expressionCol >= mt & expressionCol <= ut);
    CH = union(model.rxns(indH),core); %#ok<FNDSB>
    CM = model.rxns(indM); %#ok<FNDSB>
        cMod = MBA(model, CM, CH, tol);
        nrun = nrun + 1;
        [~,~, ~, fluxInConsistentRxnBool] = findFluxConsistentSubset(cMod,paramConsistency);
        if (sum(fluxInConsistentRxnBool) == 0)
              incon = false;
        end
    end
    if incon
        errmsg = 'Output model is inconsistent';
        save(['INC_',tName,'.mat'],'errmsg')
    end
    disp(['Number of rxns: ',num2str(numel(cMod.rxns))])
    cMod.name = tName;
    save([tName,'.mat'],tName)
    disp(' ')
end