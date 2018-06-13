function setUpOpdamSLURM(modelName, cellLine, outputDir)
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
outputRemotePrefix = '/scratch/uda2013/Opdam/';
cobraToolboxPath = '/home/uda2013/cobratoolbox/';

figures = {'U', 'C', 'S'};
methodsToGenerateModels = {'GIMME', 'fastcore', 'MBA', 'iMAT', 'INIT', 'mCADRE'};
fprintf('Setting up analysis on cell line %s\n', cellLine);

submissionFile = sprintf('allSubmitSBATCH%s_%s.sh', cellLine, modelName);
fidSubmission = fopen(submissionFile, 'w');
fprintf(fidSubmission, '#!/bin/bash\n\n');

for method = methodsToGenerateModels
    timeToRun = 48;
    if (any(strcmp(method, {'mCADRE', 'MBA', 'iMAT'})))
        for currentFig = figures
            bbs = {'B'}; % Deleted 'F' type, since it is designed (according to article) to remove the biomass function. Not sure what is the point of this.
            if(strcmp(currentFig{:}, 'C'))
                bbs = {'B', 'H'};
            end
            for currentbb = bbs
                shFileName = sprintf('subTh_%s_%s_%s_%s.sh', method{:}, modelName, cellLine, strcat(currentFig{:}, currentbb{:}));
                fid = fopen(shFileName, 'w');
                fprintf(fid, '#!/bin/bash\n#SBATCH -t 0-%d:00\n#SBATCH --ntasks=3\n#SBATCH --mem-per-cpu=2560M\n#SBATCH --output=%%x-%%j.out\n\n\n\n', timeToRun);
                fprintf(fid, 'cd %s\n', [outputRemotePrefix outputName]);
                fprintf(fid, 'matlab -nodisplay -r \"');
                fprintf(fid, 'addpath(''%s''); ', cobraToolboxPath); % So matlab on the server has cobratoolbox
                fprintf(fid, 'initCobraToolbox; ');
                fprintf(fid, 'addpath(genpath(''%s'')); ', opdamRemoteDirectory);
                fprintf(fid, 'runOpdamOnServer %s %s %s ''%s'' %s; \"\n', method{:}, currentFig{:}, currentbb{:}, modelName, cellLine);
                fclose(fid);
                fprintf(fidSubmission, 'sbatch %s\n', shFileName);
            end
        end
    else
        shFileName = sprintf('subTh_%s_%s_%s.sh', method{:}, modelName, cellLine);
        fid = fopen(shFileName, 'w');
        fprintf(fid, '#!/bin/bash\n#SBATCH -t 0-%d:00\n#SBATCH --ntasks=3\n#SBATCH --mem-per-cpu=2560M\n#SBATCH --output=%%x-%%j.out\n\n\n\n', timeToRun);
        fprintf(fid, 'cd %s\n', [outputRemotePrefix outputName]);
        fprintf(fid, 'matlab -nodisplay -r \"');
        fprintf(fid, 'addpath(''%s''); ', cobraToolboxPath); % So matlab on the server has cobratoolbox
        fprintf(fid, 'initCobraToolbox; ');
        fprintf(fid, 'addpath(genpath(''%s'')); ', opdamRemoteDirectory);
        for currentFig = figures
            bbs = {'B'};
            if (strcmp(method, 'GIMME'))
                fprintf(fid, 'runOpdamOnServer %s %s %s ''%s'' %s; ', method{:}, currentFig{:}, 'B', modelName, cellLine);
            else
                if(strcmp(currentFig{:}, 'C'))
                    bbs = {'B', 'H'};
                end
                for currentbb = bbs
                    fprintf(fid, 'runOpdamOnServer %s %s %s ''%s'' %s; ', method{:}, currentFig{:}, currentbb{:}, modelName, cellLine);
                end
            end
        end
        fprintf(fid, '"\n');
        fclose(fid);
    end
    fprintf(fidSubmission, 'sbatch %s\n', shFileName);
end

fclose(fidSubmission);

end