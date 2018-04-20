%% load model and replace HGNC ids with entrez IDs
rootDirectory = '~/Documents/MetaboGenomics';
opdamDirectory = fullfile(rootDirectory, 'Codes', 'ModelPerformanceTests', 'OpdamCellSystems2017');
recon2_2 = readCbModel('~/Documents/MetaboGenomics/Data/RECON2.2.mat');
hgncToEnsemblFile = '~/Documents/Data/General/HGNC_gene_with_protein_product.2018-Jan-22.txt';
hgncEnsembl = readtable(hgncToEnsemblFile);
hgncComplexLocusrFile = '~/Documents/Data/General/HGNCcomplex_locus_constituent.11-Feb-2018.txt';
hgncComplexLocus = readtable(hgncComplexLocusrFile);
hgncEtrez = [hgncEnsembl(:, {'hgnc_id', 'entrez_id'}); hgncComplexLocus(:, {'hgnc_id', 'entrez_id'})];
% Almost all RECON 2.2 genes have Entrez IDs, except for pseudogene(s) and
% read-through(s) which are irrelevan
[genesRemoved, indRemove] = setdiff(recon2_2.genes, hgncEtrez.hgnc_id);
recon2_2 = removeFieldEntriesForType(recon2_2, indRemove, 'genes', numel(recon2_2.genes)); 
[~, indModel, indEntrez] = intersect(recon2_2.genes, hgncEtrez.hgnc_id);
recon2_2.genes(indModel) = arrayfun(@(x) sprintf('%s', num2str(x)), hgncEtrez.entrez_id(indEntrez), 'UniformOutput', false);
recon2_2 = creategrRulesField(recon2_2);
%% constrain model
cellLines = {'HL60', 'K562', 'KBM7'}; %{'A375', 'HL60', 'K562', 'KBM7'};
for i=1:length(cellLines)
    currentCellLine = cellLines{i};
    model_c = constrainModel(recon2_2, 'C', currentCellLine);
    model_s = constrainModel(recon2_2, 'S', currentCellLine);
    model_u = constrainModel(recon2_2, 'U', currentCellLine);
    outputName = fullfile( opdamDirectory, 'Models and Data', currentCellLine, sprintf('recon2_2_%s.mat', currentCellLine));
    save(outputName, 'model_u', 'model_c', 'model_s', 'recon2_2');
end


% %% Original removeIsoGenes did not remove geneNames. This fixes that bug
% [~, genesToDelete] = setdiff(recon2_2.genes, model_u.genes);
% model_u.geneNames(genesToDelete) = '';
% [~, genesToDelete] = setdiff(recon2_2.genes, model_c.genes);
% model_c.geneNames(genesToDelete) = '';
% [~, genesToDelete] = setdiff(recon2_2.genes, model_s.genes);
% model_s.geneNames(genesToDelete) = '';
%% Run FASTCORE 
load('/Users/uridavidakavia/Documents/MetaboGenomics/Codes/ModelPerformanceTests/OpdamCellSystems2017/Models and data/A375/recon2_2_A375.mat', 'model_u', 'model_c', 'model_s', 'recon2_2');
cd(['~' filesep '' filesep 'Documents' filesep 'MetaboGenomics' filesep 'Data' filesep 'ModelPerformanceTests' filesep 'OpdamCellSystems2017' filesep 'RECON2.2'])
figures = {'U', 'C', 'S'};
cellLine = 'A375';
fprintf('Running analysis on cell line %s\n', cellLine);
for currentFig = [figures]
    runThresholdsFastCore(currentFig{:}, 'F', 'recon2_2', cellLine)
    runThresholdsFastCore(currentFig{:}, 'B', 'recon2_2', cellLine)
    if(strcmp(currentFig{:}, 'C'))
        runThresholdsFastCore(currentFig{:}, 'H', 'recon2_2', cellLine);
    end
end
%% Depletion ratio
[GenEss_score,Func_score] = Run_GenEss_Func('FastCore_recon2_2_UB1_A375', 'A375', '/Users/uridavidakavia/Documents/MetaboGenomics/Codes/ModelPerformanceTests/OpdamCellSystems2017/Models and data/A375/recon2_2_A375.mat');