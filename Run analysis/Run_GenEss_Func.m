function [mEssGenes, mNonGenes, GenEss_score,Func_score] = Run_GenEss_Func(contextSpecificModelName,CellLine, modelName)
%Run_GenEss_Func Performs analysis for Gene essentiality and Metabolic
%capabilities
%linearMOMA
%
% [GenEss_score,Func_score] = Run_GenEss_Func(model,CellLine,savefile)
%
%INPUT
% contextSpecificModelName  COBRA model structure of the extracted context-specific
%                           model
% CellLine                  Cell Line of the extracted context-specific model
%                           (A375, HL60, KBM7 or K562)
% constrainedModelFile      Filename containing the model_u constrained
%                           for this cell line
%
%OUTPUTS
% mEssP                     Copmuted Essential genes
% mNonP                     Computed Non essential genes
% GenEss_score              Computed Gene Essentiality score
% Func_score                Functionality score (pourcentage of capabilities
%                           recovered by the context-specific model)

%% Load the context-specific model to analyze
tmB = readCbModel(contextSpecificModelName);

%% load gene essentiality from in vitro screen
global  depl_p

if (nargin < 3)
    load(['model_u_',CellLine,'.mat']);
    load(['gene_expr_u_',CellLine,'.mat']);
    load(['gene_id_u_',CellLine,'.mat']);
    depl_p = findModelDepletionGenes_loc(model_u,['depletion_ratios_',CellLine,'.xls']);
else
    depl_p = findModelDepletionGenes_loc(tmB,['depletion_ratios_',CellLine,'.xls']);
end
%% Check gene-essentiality
% Growth rate for essential genes after KO is at most 0.01 of WT -
% TODO should rate be different?
maxgr=0.01;
tmBcRed = tmB;
if (all(tmBcRed.c == 0))
    warning('No objective reaction %s was found in model %s!\n', tmBcRed.biomassRxnAbbr);
    mEssGenes = NaN;
    mNonGenes = NaN;
    GenEss_score = NaN;
    Func_score = NaN;
    return;
end
%TODO - think about if lower bound for biomass should be zero
tmBcRed = changeRxnBounds(tmBcRed,tmBcRed.rxns(tmBcRed.c == 1),0,'l');
grRatios = singleGeneDeletion(tmBcRed,'FBA',depl_p.genes);

essentialGemeInd = grRatios <= maxgr;
nonEssentialGemeInd = grRatios > maxgr;
mEssP = depl_p.values(essentialGemeInd);
mNonP = depl_p.values(nonEssentialGemeInd);

mEssGenes = depl_p.genes(essentialGemeInd);
mNonGenes = depl_p.genes(nonEssentialGemeInd);

try
GenEss_score = ranksum(mEssP,mNonP,'tail','left');
catch ME
warning('ranksum failed with error');
warning(ME.message)
GenEss_score = NaN;
end
% Not really clear what testFunctionality is trying to do. Should
% probably read original paper
[totTest, numTests] = testFunctionality(tmB);
Func_score=totTest/numTests*100;
end

function [totTest, numTests] = testFunctionality(model)
    [~,bmMets,~] = xlsread('biomass_rxn.xls');
    bmMets = setdiff(bmMets,{'atp[c]','adp[c]','pi[c]','h2o[c]','h[c]'});
    numTests = numel(bmMets) + 1;
    if (isfield(model, 'comps'))
        comp = strcat('[', model.comps, ']');
    else
        comp = {'[c]','[e]','[g]','[l]','[m]','[n]','[r]','[x]'}; % if model doesn't have compartments
    end
    totTest = 0;
    for i = 1:numel(bmMets)
        cMet = bmMets{i}(1:end-3);
        cTest = 0;
        for j = 1:numel(comp)
            met = [cMet,comp{j}];
            if ~isempty(intersect(model.mets,met))
                nm = addReaction(model,['BMS_',met],{met},-0.5,false,0,1000,0,'Biomass Metabolite Sink');
                nm = changeObjective(nm, ['BMS_',met]);
                sol = optimizeCbModel(nm);
                if sol.f > 1e-8
                    cTest = 1;
                    break;
                end
            end
        end
        totTest = cTest + totTest;
    end
    
    warning off all
    nm = addReaction(model,'BMS_atp(c)',{'atp[c]','h2o[c]','adp[c]','pi[c]','h[c]'},[-0.5,-0.5,0.5,0.5,0.5],false,0,1000,0,'Biomass Metabolite Sink'); %Make a copy of ATP demand
    atpDMInd = strncmp(model.rxns, 'DM_atp', 6) | strcmp(model.rxns, 'ATPM'); 
    nm = changeRxnBounds(nm, model.rxns(atpDMInd) ,0,'l');
    warning on all
    nm = changeObjective(nm, 'BMS_atp(c)');
    sol = optimizeCbModel(nm);
    if sol.f > 1e-8
        totTest = totTest + 1;
    end
    
    disp(['Functionality test: ',num2str(totTest),'/',num2str(numel(bmMets)+1)])
end

%%
function genes_ratios = findModelDepletionGenes_loc(model,filename)

    % Load depletion (CRISPR) and ID data
    [~,~,raw] = xlsread(filename);
    fid = fopen('gecko_id_to_entrez.txt', 'r'); % TODO - THis file was done in XLS, and screwed up several gene names. Can probably be replaced by processing gene_info, or some newer version from the GECKO website.
    % TODO - use the CERES correction of Gecko, but not at first
    % https://figshare.com/articles/CERES_-_Meyers_Bryan_et_al_Nature_Genetics_2017_/5319388
    foo = textscan(fid, '%s %s');
    fclose(fid);
    gGeck = foo{1};
    gEntr = foo{2};
    clear foo
    gDepl = raw(:,1);
    for i = 1:length(gDepl)
        if isnumeric(gDepl{i})
            gDepl{i} = num2str(gDepl{i});
        end
    end
   
    dRatio = raw(:,2);    
    % Only get the genes which can be translated from gecko to entrez
    [~, indEntrez, indDepl] = intersect(gGeck, gDepl);
    gEntr = gEntr(indEntrez);
    dRatio = cell2mat(dRatio(indDepl));
    
    % Only get the genes which we have depletion data for and are in the
    % model
    [members, geneInds] = ismember(model.genes, gEntr);
    mRatio = dRatio(geneInds(geneInds~=0));
    mGenes = model.genes(members);

    genes_ratios.values = mRatio;
    genes_ratios.genes = mGenes;
    
end