function runThresholdsGIMME(figName, modelName, cellLine, overWrite)

if (nargin < 4 || isempty(overWrite))
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

figName = [figName,'B'];
singleRun(model, expressionCol, ths.p10, 0.9, figName, 1, modelName, cellLine, tol, overWrite)
singleRun(model, expressionCol, ths.mean, 0.9, figName, 2, modelName, cellLine, tol, overWrite)
singleRun(model, expressionCol, ths.p25, 0.9, figName, 3, modelName, cellLine, tol, overWrite)
singleRun(model, expressionCol, ths.p50, 0.9, figName, 4, modelName, cellLine, tol, overWrite)
end

function singleRun(model, expressionCol, ut, obj_frac, figName, id, modelName, cellLine, tol, overWrite)
    tName = ['GIMME_', cellLine, '_', figName, num2str(id), '_', modelName];
    disp(tName)
    disp('RUNNING GIMME...')
    optionsLocal = struct('solver', 'GIMME', 'expressionRxns', {expressionCol},...
        'threshold', ut, 'obj_frac', obj_frac);
    paramsLocal.epsilon = tol;
    if (~overWrite && exist([pwd '/' tName '.mat'], 'file') ~=2)
        overWrite = 1;
    end
    if (overWrite)
        try
            cMod = createTissueSpecificModel(model, optionsLocal, 1, [], paramsLocal);
            cMod.name = tName;
            writeCbModel(cMod, 'mat', tName);
        catch ME
            warning('Failed to run GIMME on model %s, figure %s with cell line %s', modelName, [figName num2str(id)], cellLine);
            warning(ME.message)
        end
    end
end