function modelNew=removeIsoGenes(model)
    % Gets rid of the Entrez suffices ".x" in the gene identifiers in the
    % entire model, also updates GPR rules.
    % Input: model with genes in the form "entrezID.x"
    % Output: model where all genes are in the form "enrezID"
    
    model = rmfield(model, {'rules', 'rxnGeneMat'});
    model.genes=[];
    model.grRules = regexprep(model.grRules, '[.]\d*', '');
    warning off all
    model = buildRxnGeneMat(model);
    warning on all
    modelNew=model;
end