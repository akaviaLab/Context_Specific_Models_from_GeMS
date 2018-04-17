function setUpOpdamMOAB(modelName, cellLine, outputDir)
% This sets up a run of Opdam 2017 (making context-sepcific models using
% GIMME, FASTCORE, MBA, iMAT, INIT, MBA and mCadre) on a specific cell line
% (out of 4 options: 'A375', 'HL60', 'K562', 'KBM7') given as the second
% input.
% The run is set up on the MOAB scheduler, and uses Torque commands

if (nargin < 2)
    fprintf("Usage: setupOpdamMOAB(model, cellLine, outputDir)\n");
    fprintf(" .     setupOPdamMOAB([],    'A375',   'output')\n");
    return;
end

cd(outputDir);
if (isempty(modelName))
    outputName = 'Recon1';
else
    outputName = modelName;
end
mkdir(outputName);
cd(outputName);

opdamRemoteDirectory = '/home/uda2013/OpdamCellSystems2017';

figures = {'U', 'C', 'S'};
methodsToGenerateModels = {'GIMME', 'fastcore', 'MBA', 'iMAT', 'INIT', 'mCADRE'};
fprintf('Setting up analysis on cell line %s\n', cellLine);
for method = methodsToGenerateModels
    for currentFig = figures
        bbs = {'F', 'B'};
        if(strcmp(currentFig{:}, 'C'))
            bbs = [bbs, 'H'];
        end
        for currentbb = bbs
            shFileName = sprintf('subThresholds%s_%s_%s%s.sh', method{:}, modelName, cellLine, strcat(currentFig{:}, currentbb{:}));
            %matlabFileName = sprintf('runThresholds%s_%s_%s%s.m', method{:}, modelName, cellLine, strcat(currentFig{:}, currentbb{:}));
            fid = fopen(shFileName, 'w');
            fprintf(fid, '#PBS -l walltime=12:00:00\n#PBS -l nodes=1:ppn=12\n#PBS -A iii-973-aa\n\n\n');
            fprintf(fid, 'cd /gs/scratch/uda2013/Opdam/%s\n', outputName);
            fprintf(fid, 'matlab -nodisplay -r \"');
            fprintf(fid, 'addpath(genpath(''%s''); ', opdamRemoteDirectory);
            fprintf(fid, 'runOpdamOnServer %s %s %s ''%s'' %s;', method{:}, currentFig{:}, currentbb{:}, modelName, cellLine);            
            fclose(fid);
        end
    end
end

end