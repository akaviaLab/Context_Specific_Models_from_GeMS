function runOpdamOnServer(extractionMethod, figName, bbLetter, modelName, cellLine, overwrite)
changeCobraSolver('ibm_cplex', 'all') % Glpk fails when using CheckModelConsistency (on the server?)
% give number of workers for parallelization
% pc = parcluster('local');
% mkdir(strcat('/localscratch/', getenv('SLURM_JOBID')));
% pc.JobStorageLocation = strcat('/localscratch/', getenv('SLURM_JOBID')); % Which directory are the threads using to communicate
% nworkers = str2double(getenv('SLURM_JOB_NUM_NODES')); %PBS_NUM_PPN
% parpool(pc, nworkers) % starts the parallel pool

if nargin < 6
    overwrite=1;
end

numberModels = 4;

switch extractionMethod
    case 'GIMME'
        runThresholdsGIMME(figName, modelName, cellLine, overwrite);
    case 'fastcore'
        runThresholdsFastCore(figName, bbLetter, modelName, cellLine, overwrite);
    case 'iMAT'
        runThresholdsiMAT(figName, bbLetter, modelName, cellLine, overwrite);
        numberModels = 10;
    case 'INIT'
        runThresholdsINIT(figName, bbLetter, modelName, cellLine, overwrite);
    case 'MBA'
        runThresholdsMBA(figName, bbLetter, modelName, cellLine, overwrite);
        numberModels = 10;
    case 'mCADRE'
        runThresholdsmCADRE(figName, bbLetter, modelName, cellLine, overwrite);
end

end
