%% load model and replace HGNC ids with entrez IDs
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
model_c = constrainModel(recon2_2, 'C', 'A375');
model_s = constrainModel(recon2_2, 'S', 'A375');
model_u = constrainModel(recon2_2, 'U', 'A375');

save('/Users/uridavidakavia/Documents/MetaboGenomics/Codes/ModelPerformanceTests/OpdamCellSystems2017/Models and data/A375/recon2_2_A375.mat', 'model_u', 'model_c', 'model_s', 'recon2_2');
% %% Original removeIsoGenes did not remove geneNames. This fixes that bug
% [~, genesToDelete] = setdiff(recon2_2.genes, model_u.genes);
% model_u.geneNames(genesToDelete) = '';
% [~, genesToDelete] = setdiff(recon2_2.genes, model_c.genes);
% model_c.geneNames(genesToDelete) = '';
% [~, genesToDelete] = setdiff(recon2_2.genes, model_s.genes);
% model_s.geneNames(genesToDelete) = '';
%% Run FASTCORE 
load('/Users/uridavidakavia/Documents/MetaboGenomics/Codes/ModelPerformanceTests/OpdamCellSystems2017/Models and data/A375/recon2_2_A375.mat', 'model_u', 'model_c', 'model_s');
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
