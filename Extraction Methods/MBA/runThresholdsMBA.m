function runThresholdsMBA(figName, bb, modelName, cellLine, overWrite)

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
    %UNCONSTRAINED
    tol = 1e-6;
    expressionCol = mapExpressionToReactions(model_u, expressionData_u);
    model = model_u;
end

if strcmp(figName,'C')
    %CONSTRAINED
    tol = 1e-8;
    expressionCol = mapExpressionToReactions(model_c, expressionData_c);
    model = model_c;
end
if strcmp(figName,'S')
    %SEMI-CONSTRAINED
    tol = 1e-6;
    expressionCol = mapExpressionToReactions(model_s, expressionData_s);
    model = model_s;
end

core = {};
if strcmp(bb, 'B')
    biomassRxnInd = strncmpi(model.rxns, 'biomass', 7);
    biomassRxn = model.rxns(biomassRxnInd);
    atpDMInd = strncmp(model.rxns, 'DM_atp', 6) | strcmp(model.rxns, 'ATPM');
    model = changeRxnBounds(model, biomassRxn, blb, 'l'); %Force biomass and ATP demand to be active
    core = model.rxns(biomassRxnInd | atpDMInd);
    figName = [figName,'B'];
end
if strcmp(bb, 'F')
    model = addBiomassSinks(model);
    figName = [figName,'F'];
end
if strcmp(bb, 'H')
    biomassRxnInd = strncmpi(model.rxns, 'biomass', 7);
    biomassRxn = model.rxns(biomassRxnInd);
    atpDMInd = strncmp(model.rxns, 'DM_atp', 6) | strcmp(model.rxns, 'ATPM');
    model = changeRxnBounds(model, biomassRxn, 1e-3, 'l'); %Force biomass and ATP demand to be active
    core = model.rxns(biomassRxnInd | atpDMInd);
    figName = [figName,'H'];
end
run_MBA(core, model, expressionCol, figName, ths.p10, ths.p10, 1, modelName, tol, cellLine, overWrite)
run_MBA(core, model, expressionCol, figName, ths.mean, ths.mean, 2, modelName, tol, cellLine, overWrite)
run_MBA(core, model, expressionCol, figName, ths.p25, ths.p25, 3, modelName, tol, cellLine, overWrite)
run_MBA(core, model, expressionCol, figName, ths.p50, ths.p50, 4, modelName, tol, cellLine, overWrite)
run_MBA(core, model, expressionCol, figName, ths.mean, ths.p10, 5, modelName, tol, cellLine, overWrite)
run_MBA(core, model, expressionCol, figName, ths.p25, ths.p10, 6, modelName, tol, cellLine, overWrite)
run_MBA(core, model, expressionCol, figName, ths.p50, ths.p10, 7, modelName, tol, cellLine, overWrite)
run_MBA(core, model, expressionCol, figName, ths.p25, ths.mean, 8, modelName, tol, cellLine, overWrite)
run_MBA(core, model, expressionCol, figName, ths.p50, ths.mean, 9, modelName, tol, cellLine, overWrite)
run_MBA(core, model, expressionCol, figName, ths.p50, ths.p25, 10, modelName, tol, cellLine, overWrite)
end

function run_MBA(core, model, expressionCol, figName, mt, ut, id, modelName, tol, cellLine, overWrite)
    tName = ['MBA', '_', figName, num2str(id), '_', cellLine, '_', modelName];
    disp(tName)
    
    indH = find(expressionCol > ut);
    indM = find(expressionCol >= mt & expressionCol <= ut);
    CH = union(model.rxns(indH),core); %#ok<FNDSB>
    CM = model.rxns(indM); %#ok<FNDSB>
    NC = setdiff(model.rxns, union(CH, CM));
    %Biomass metabolite sinks are not pruned in call_MBA
    biomassMetSinksInd = contains(NC, 'BMS_'); 
    CM = union(CM, NC(biomassMetSinksInd));
    
    if (~overWrite && exist([pwd '/' tName '.mat'], 'file') ~=2)
        overWrite = 1;
    end
    if (overWrite)
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
        
end