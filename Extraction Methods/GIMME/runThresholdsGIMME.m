function runThresholdsGIMME(figName, modName)
    %initCobraToolbox
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
        tol = 1e-6;
        figName = [figName,'B'];
        disp('UNCONSTRAINED MODEL')
        expressionCol = mapExpressionToReactions(model_u, expressionData_u);
        singleRun(model_u, expressionCol, ths.p10, 0.9, figName, 1, modName, tol)
        singleRun(model_u, expressionCol, ths.mean, 0.9, figName, 2, modName, tol)
        singleRun(model_u, expressionCol, ths.p25, 0.9, figName, 3, modName, tol)
        singleRun(model_u, expressionCol, ths.p50, 0.9, figName, 4, modName, tol)
    end
    if strcmp(figName,'C')
        tol = 1e-8;
        disp('CONSTRAINED MODEL')
        figName = [figName,'B'];
        expressionCol = mapExpressionToReactions(model_c, expressionData_c);
        singleRun(model_c, expressionCol, ths.p10, 0.9, figName, 1, modName, tol)
        singleRun(model_c, expressionCol, ths.mean, 0.9, figName, 2, modName, tol)
        singleRun(model_c, expressionCol, ths.p25, 0.9, figName, 3, modName, tol)
        singleRun(model_c, expressionCol, ths.p50, 0.9, figName, 4, modName, tol)
    end
    if strcmp(figName,'S')
        tol = 1e-6;
        disp('SEMI-CONSTRAINED MODEL')
        figName = [figName,'B'];
        expressionCol = mapExpressionToReactions(model_s, expressionData_s);
        singleRun(model_s, expressionCol, ths.p10, 0.9, figName, 1, modName, tol)
        singleRun(model_s, expressionCol, ths.mean, 0.9, figName, 2, modName, tol)
        singleRun(model_s, expressionCol, ths.p25, 0.9, figName, 3, modName, tol)
        singleRun(model_s, expressionCol, ths.p50, 0.9, figName, 4, modName, tol)
    end
end

function singleRun(model, expressionCol, ut, obj_frac, figName, id, modName, tol)
    tName = ['GIMME_',figName,num2str(id),'_',modName];
    disp(tName)
    disp('RUNNING GIMME...')
    cMod= call_GIMME(model, expressionCol, ut, obj_frac, tol);
    cMod.name = tName;
    eval([tName,'= addMinGeneField(cMod,gene_exp,gene_names,ut,0)']);
    save([tName,'.mat'],tName)
end