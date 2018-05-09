function runThresholdsGIMME(figName, modelName, cellLine)
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
        figName = [figName,'B'];
        disp('UNCONSTRAINED MODEL')
        expressionCol_u = mapExpressionToReactions(model_u, expressionData_u);
        singleRun(model_u, expressionCol_u, ths.p10, 0.9, figName, 1, modelName, cellLine, tol)
        singleRun(model_u, expressionCol_u, ths.mean, 0.9, figName, 2, modelName, cellLine, tol)
        singleRun(model_u, expressionCol_u, ths.p25, 0.9, figName, 3, modelName, cellLine, tol)
        singleRun(model_u, expressionCol_u, ths.p50, 0.9, figName, 4, modelName, cellLine, tol)
    end
    if strcmp(figName,'C')
        tol = 1e-8;
        disp('CONSTRAINED MODEL')
        figName = [figName,'B'];
        expressionCol_c = mapExpressionToReactions(model_c, expressionData_c);
        singleRun(model_c, expressionCol_c, ths.p10, 0.9, figName, 1, modelName, cellLine, tol)
        singleRun(model_c, expressionCol_c, ths.mean, 0.9, figName, 2, modelName, cellLine, tol)
        singleRun(model_c, expressionCol_c, ths.p25, 0.9, figName, 3, modelName, cellLine, tol)
        singleRun(model_c, expressionCol_c, ths.p50, 0.9, figName, 4, modelName, cellLine, tol)
    end
    if strcmp(figName,'S')
        tol = 1e-6;
        disp('SEMI-CONSTRAINED MODEL')
        figName = [figName,'B'];
        expressionCol_s = mapExpressionToReactions(model_s, expressionData_s);
        singleRun(model_s, expressionCol_s, ths.p10, 0.9, figName, 1, modelName, cellLine, tol)
        singleRun(model_s, expressionCol_s, ths.mean, 0.9, figName, 2, modelName, cellLine, tol)
        singleRun(model_s, expressionCol_s, ths.p25, 0.9, figName, 3, modelName, cellLine, tol)
        singleRun(model_s, expressionCol_s, ths.p50, 0.9, figName, 4, modelName, cellLine, tol)
    end
end

function singleRun(model, expressionCol, ut, obj_frac, figName, id, modelName, cellLine, tol)
    tName = ['GIMME_',modelName, '_', figName, num2str(id),'_',cellLine];
    disp(tName)
    disp('RUNNING GIMME...')
    cMod= call_GIMME(model, expressionCol, ut, obj_frac, tol);
    writeCbModel(cMod, 'mat', tName);
    optionsLocal = struct('solver', 'GIMME', 'expressionRxns', {expressionCol},...
        'threshold', ut, 'obj_frac', obj_frac);
    paramsLocal.epsilon = tol;
    try
        cMod2 = createTissueSpecificModel(model, optionsLocal, 1, [], paramsLocal);
        cMod2.name = tName;
        writeCbModel(cMod2, 'mat', [tName '_2']);
        if (~isSameCobraModel(cMod, cMod2))
            warning('In GIMME model %s, cell line %s, fig %s, id %d, the old and new models are different!\n', modelName, cellLine, figName, id);
        end
    catch ME
        warning('Failed to run GIMME on model %s, figure %s with cell line %s', modelName, [figName num2str(id)], cellLine);
        warning(ME.message)
    end
end