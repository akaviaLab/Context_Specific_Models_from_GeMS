function modelNew=removeIsoGenes(model)
    % Gets rid of the Entrez suffices ".x" in the gene identifiers in the
    % entire model, also updates GPR rules.
    % Input: model with genes in the form "entrezID.x"
    % Output: model where all genes are in the form "enrezID"
    
    model = rmfield(model, {'rules', 'rxnGeneMat'});
    relevantFields = getModelFieldsForType(model, 'genes');
    relevantFields = setdiff(relevantFields, 'genes');
    originalGenes = model.genes;
    originalGenes = regexprep(originalGenes, '[.]\d*', '');
    model.genes = [];
    model.grRules = regexprep(model.grRules, '[.]\d*', '');
    warning off all
    model = buildRxnGeneMat(model);
    warning on all
    modelNew=model;
    [~, indToKeep, indNew] = intersect(originalGenes, modelNew.genes);
    indToRemove = setdiff(1:length(originalGenes), indNew);
    for i=1:length(relevantFields)
        modelNew.(relevantFields{i})(indNew) = model.(relevantFields{i})(indToKeep);
        modelNew.(relevantFields{i})(indToRemove) = '';
    end
end