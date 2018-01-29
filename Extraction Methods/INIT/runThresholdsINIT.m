function runThresholdsINIT(figName, bb, modName)
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
    runtime = 3600;
    if strcmp(bb,'B')
        model = changeRxnBounds(model_u, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
        figName = [figName,'B'];
    end
    if strcmp(bb,'F')
        model = addBiomassSinks(model_u);
        bmsInd=~cellfun(@isempty,strfind(model.rxns,'BMS_'));
        figName = [figName,'F'];
    end
    expressionCol = mapExpressionToReactions(model_u, expressionData_u);
end

if strcmp(figName,'C')
    tol = 1e-8;
    runtime = 7200;
    if strcmp(bb,'B')
        model = changeRxnBounds(model_c, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
        figName = [figName,'B'];
    end
    if strcmp(bb,'F')
        model = addBiomassSinks(model_c);
        bmsInd=~cellfun(@isempty,strfind(model.rxns,'BMS_'));
        figName = [figName,'F'];
    end
    if strcmp(bb,'H')
        model = changeRxnBounds(model_c, 'Biomass_reaction', 1e-3, 'l'); %Force biomass and ATP demand to be active
        figName = [figName,'H'];
    end
    expressionCol = mapExpressionToReactions(model_c, expressionData_c);
end

if strcmp(figName,'S')
    tol = 1e-6;
    runtime = 3600;
    if strcmp(bb,'B')
        model = changeRxnBounds(model_s, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
        figName = [figName,'B'];
    end
    if strcmp(bb,'F')
        model = addBiomassSinks(model_s);
        bmsInd=~cellfun(@isempty,strfind(model.rxns,'BMS_'));
        figName = [figName,'F'];
    end
    expressionCol = mapExpressionToReactions(model_s, expressionData_s);
end


%w1
w1 = zeros(length(expressionCol),1);
w1(expressionCol >= 0) = 5*log(expressionCol(expressionCol>=0)/ths.p10);
w1(expressionCol < 0) = -2; % "unknown" entries get a weight of -2
w1(w1 < -max(w1)) = -max(w1);
if strcmp(bb,'B') || strcmp(bb,'H')
    w1(ismember(model.rxns, {'Biomass_reaction','DM_atp(c)'})) = max(w1); %Biomass and ATP demand get high weight
else
    w1(bmsInd) = 0;
end

%w2
w2 = zeros(length(expressionCol),1);
w2(expressionCol >= 0) = 5*log(expressionCol(expressionCol>=0)/ths.mean);
w2(expressionCol < 0) = -2; % "unknown" entries get a weight of -2
w2(w2 < -max(w2)) = -max(w2);
if strcmp(bb,'B') || strcmp(bb,'H')
    w2(ismember(model.rxns, {'Biomass_reaction','DM_atp(c)'})) = max(w2); %Biomass and ATP demand get high weight
end

%w3
w3 = zeros(length(expressionCol),1);
w3(expressionCol >= 0) = 5*log(expressionCol(expressionCol>=0)/ths.p25);
w3(expressionCol < 0) = -2; % "unknown" entries get a weight of -2
w3(w3 < -max(w3)) = -max(w3);
if strcmp(bb,'B') || strcmp(bb,'H')
    w3(ismember(model.rxns, {'Biomass_reaction','DM_atp(c)'})) = max(w3); %Biomass and ATP demand get high weight
else
    w3(bmsInd) = 0;
end

%w4
w4 = zeros(length(expressionCol),1);
w4(expressionCol >= 0) = 5*log(expressionCol(expressionCol>=0)/ths.p50);
w4(expressionCol < 0) = -2; % "unknown" entries get a weight of -2
w4(w4 < -max(w4)) = -max(w4);
if strcmp(bb,'B') || strcmp(bb,'H')
    w4(ismember(model.rxns, {'Biomass_reaction','DM_atp(c)'})) = max(w4); %Biomass and ATP demand get high weight
else
    w4(bmsInd) = 0;
end
    
if strcmp(figName,'UB') || strcmp(figName,'UF')
    %UNCONSTRAINED
    epsil = 1;
    run_INIT(model,ths.p25, figName, epsil, w1, 1, modName, tol, runtime);
    run_INIT(model,ths.p50, figName, epsil, w2, 2, modName, tol, runtime);
    run_INIT(model,ths.p10, figName, epsil, w3, 3, modName, tol, runtime);
    run_INIT(model,ths.mean, figName, epsil, w4, 4, modName, tol, runtime);
end

if strcmp(figName,'CB') || strcmp(figName,'CF') || strcmp(figName,'CH')
    %CONSTRAINED
    epsil = 1e-6;
    run_INIT(model,ths.p25, figName, epsil, w1, 1, modName, tol, runtime);
    run_INIT(model,ths.p50, figName, epsil, w2, 2, modName, tol, runtime);
    run_INIT(model,ths.p10, figName, epsil, w3, 3, modName, tol, runtime);
    run_INIT(model,ths.mean, figName, epsil, w4, 4, modName, tol, runtime);
end
if strcmp(figName,'SB') || strcmp(figName,'SF')
    %SEMI-CONSTRAINED
    epsil = 1;
    run_INIT(model,ths.p25, figName, epsil, w1, 1, modName, tol, runtime);
    run_INIT(model,ths.p50, figName, epsil, w2, 2, modName, tol, runtime);
    run_INIT(model,ths.p10, figName, epsil, w3, 3, modName, tol, runtime);
    run_INIT(model,ths.mean, figName, epsil, w4, 4, modName, tol, runtime);

end
exit;
end

function run_INIT(model, ~, figName, ~, w, id, modName, tol, runtime)
    tName = ['INIT_',figName, num2str(id),'_',modName];
    disp(tName)
    optionsLocal = struct('weights', w, 'tol', tol, 'runtime', runtime');
    optionsLocal.solver = 'INIT';
    optionsLocal.logfile = [tName,'.txt'];
    cMod = createTissueSpecificModel(model, optionsLocal, 1);
    cMod.name = tName;
    save([tName,'.mat'],tName)
end