% This script will compare the results between the published code, and the
% code using the COBRATOOLBOX

% Summary of results for A375
% GIMME GIMME__CB3_A375 seems to be different in size of A. Why? Because
% the published version uses an old version of convertToIrreversible. If
% using the toolbox version, it works.
% FASTCORE - ??
% INIT & iMAT - the default epsilon in createTissueModel works for U and S,
% not for C. Need to consider adding epsilon as an option.
% MBA is massively different, need to check why - fastcc in Opdam is old
% and probably wrong
% mCADRE in Opdam is weird, since it has hardcoded checkprecursor
% metaoblites, using the toolbox for all future work

methodsToCheck = {'GIMME', 'iMAT', 'INIT', 'MBA', 'FastCore'};
figures = {'C', 'U', 'S'};

fprintf('\tFigure\tBB\t#checked\t#match\n');
for i=1:length(methodsToCheck)
    curMethod = methodsToCheck{i};
    for j=1:length(figures)
        curFigure = figures{j};
        bbs = {'B', 'F'};
        if (strcmp(curFigure, 'C'))
            bbs = [bbs, 'H'];
        end
        if (strcmp(curMethod, 'GIMME'))
            bbs = {'B'};
        end
        for k=1:length(bbs)
            curBB = bbs{k};
            fprintf('%s\t%s\t%s\t', curMethod, curFigure, curBB);
            currentFileList = dir([curMethod, '*', curFigure '*', curBB, '*.mat']);
            fprintf('%d\t', length(currentFileList));
            numModelsMatching = 0;
            indFile = 1;
            mismatchingFields = cell(0);
            while (indFile < length(currentFileList))
                load(currentFileList(indFile).name);
                model1 = model;
                load(currentFileList(indFile+1).name);
                model2 = model;
                indFile = indFile + 2;
                try
                    [a, b, c] = isSameCobraModel(model1, model2);
                    if (a)
                        numModelsMatching = numModelsMatching + 1;
                    else
                        mismatchingFields = union(mismatchingFields, c(logical(b)));
                    end
                catch
                end
            end
            fprintf('%d\t', numModelsMatching * 2);
            fprintf('%s,', mismatchingFields{:});
            fprintf('\n');
        end
    end
end