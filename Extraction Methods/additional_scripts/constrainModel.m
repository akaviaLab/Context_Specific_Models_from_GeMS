function new_model = constrainModel(model, constr,cellType)
    % Constrains and modifies Recon1 using the exometabolomic data
    % Input:
    %   constr - type of constraint that is used on exhcanges: "U", "S" or "C"
    %         U: Unconstrained, only constrain carbon-source metabolites to 
    %            a maximum uptake of 10
    %         S: Semi-constrained, constrain the direction of measured
    %            exchange reaction and do not allow uptake of other
    %            carbon-source metabolites
    %         C: Consrained, constrain using the exometabolomic values and 
    %            do not allow uptake of other carbon-source metabolites
    % cellType - name of the cell type: "A375" or "HL60"
    % Output:
    %   new_model - constrained and modified version of recon 1. This model
    %               is consistent and can be used as input model for the
    %               MEMs. Units are in mmol/gDw/h
    %
    % This script expects the following files to be present:
    %  - exometabolomics_minmax_<cellType>.xlsx
    
    new_model = findSExRxnInd(model);
    
    if strcmp(constr,'C')
        tol = 1e-8;
    elseif strcmp(constr,'U') || strcmp(constr,'S')
        tol = 1e-6;
    else
        error('Unknown type of constraint, use "U", "S" or "C"')
    end
    
    %Initially constrain exchange reactions to zero uptake for S and C
    if strcmp(constr,'S') || strcmp(constr,'C')
        rxnsOpen = {'EX_ca2(e)';'EX_cl(e)';'EX_co(e)';'EX_co2(e)';'EX_fe2(e)';'EX_fe3(e)';'EX_h(e)';'EX_h2o(e)';'EX_h2o2(e)';'EX_hco3(e)';'EX_i(e)';'EX_k(e)';'EX_na1(e)';'EX_nh4(e)';'EX_no(e)';'EX_o2(e)';'EX_o2s(e)';'EX_oh1';'EX_oxa(e)';'EX_pi(e)';'EX_sel(e)';'EX_so4(e)';'EX_tcynt(e)';'EX_tsul(e)'};
        new_model = changeRxnBounds(new_model,new_model.rxns(new_model.ExchRxnBool),0,'l');
        new_model = changeRxnBounds(new_model,rxnsOpen,-1000,'l');
    end
    
    %Now set exchange reaction rates from measured values
    [exVal,exName,~] = xlsread(['exometabolomics_minmax_',cellType,'.xlsx']);
    essRxns = {'EX_his_L(e)';'EX_ocdca(e)';'EX_ocdcea(e)';'EX_Tyr_ggn(e)';'EX_cys_L(e)'}; %Reactions not in dataset that are allowed uptake
    %Constrain quantitatively
    if strcmp(constr,'C')
        new_model = changeRxnBounds(new_model,exName,exVal(:,1),'l');
        new_model = changeRxnBounds(new_model,exName,exVal(:,2),'u');
        new_model = changeRxnBounds(new_model,essRxns,-1000*ones(numel(essRxns),1),'l');
        [~,maxf] = fluxVariability(new_model,100,'max',essRxns);
        new_model = changeRxnBounds(new_model,essRxns,maxf,'l');
    end
    %Constrain qualitatively (direction)
    if strcmp(constr,'S')
        minf = exVal(:,1);
        maxf = exVal(:,2);
        minf(minf<0) = -1000;
        maxf(maxf>0) = 1000;
        new_model = changeRxnBounds(new_model,exName,minf,'l');
        new_model = changeRxnBounds(new_model,exName,maxf,'u');
        new_model = changeRxnBounds(new_model,essRxns,-1000*ones(numel(essRxns),1),'l');
    end
    %Constrain carbon-source metabolites for U and S
    if strcmp(constr,'U') || strcmp(constr,'S')
        orgExRxns = find_organic_ex_rxns(new_model);
        for i = 1:numel(orgExRxns)
            ori = find(ismember(new_model.rxns, orgExRxns{i}));
            if new_model.lb(ori) < -10
                new_model.lb(ori) = -10;
            end
        end
    end
    
    %Remove inconsistent reactions
    %activeRxns = fastcc(new_model, tol);
    %inactiveRxns = new_model.rxns(setdiff(1:length(new_model.rxns), activeRxns));
    %new_model = removeRxns(new_model, inactiveRxns);
    
    %Remove isoforms to make genes unique and also remove genes no longer
    %mapping to any reaction
    new_model = removeIsoGenes(new_model);
    
    %Check biomass flux
    result = optimizeCbModel(new_model);
    disp(['Maximum biomass rate: ',num2str(result.f)])
    if result.f < tol
       warning('Biomass as objective cannot carry flux')
    end
end

function organicExRxns = find_organic_ex_rxns(model)
    %Sub-function taken from mCADRE to identify carbon-source metabolites

    ind2EX = ~cellfun(@isempty,strfind(model.rxns,'EX_'));
    exRxns = model.rxns(ind2EX);

    % Note: should add a warning if the metFormulas field is empty
    % Organic metabolites are defined as those containing carbon (C) and hydrogen
    % (H); these are identified by checking molecular formulas
    if isempty(model.metFormulas); warning('metFormulas field is empty'); end
    carbonMets = ~cellfun('isempty', regexp(model.metFormulas, 'C'));
    hydrogenMets = ~cellfun('isempty', regexp(model.metFormulas, 'H'));
    is_organic = carbonMets & hydrogenMets;
    organicMets = model.mets(is_organic);

    organicRxns = findRxnsFromMets(model, organicMets);
    organicExRxns = intersect(organicRxns, exRxns);

    % The following reactions exchange organic metabolites (e.g., R-groups that
    % comprise lipid tails), but don't contain H in their specified formulas; OR,
    % as in the case of Tyr-ggn, include protein compounds
    organicExRxns = [organicExRxns; ...
        'EX_Rtotal(e)'; 'EX_Rtotal2(e)'; 'EX_Rtotal3(e)'; 'EX_Tyr_ggn(e)';'EX_peplys(e)']; % ; ...
        % 'UP_Tyr_ggn[c]']; % This rxn is not in Recon 1 - may be something specific
        % to Recon 2...
    organicExRxns = setdiff(organicExRxns, 'EX_hco3(e)');
end