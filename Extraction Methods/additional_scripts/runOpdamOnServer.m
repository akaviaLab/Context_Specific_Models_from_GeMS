function runOpdamOnServer(extractionMethod, figName, bbLetter, modelName, cellLine)

addpath('/gs/home/uda2013/cobratoolbox/'); % So matlab on the server has cobratoolbox
% give number of workers for parallelization
pc = parcluster('local');
pc.JobStorageLocation = strcat('/localscratch/', getenv('PBS_JOBID')); % Which directory are the threads using to communicate
nworkers = str2double(getenv('PBS_NUM_PPN'));
parpool(pc, nworkers) % starts the parallel pool

switch extractionMethod
    case 'iMAT'
        runThresholdsiMAT(figName, bbLetter, modelName, cellLine);
    case 'GIMME'
        runThresholdsGIMME(figName, modelName, cellLine);
    case 'INIT'
        runThresholdsINIT(figName, bbLetter, modelName, cellLine);
    case 'MBA'
        runThresholdsMBA(figName, bbLetter, modelName, cellLine);
    case 'mCADRE'
        runThresholdsmCADRE(figName, bbLetter, modelName, cellLine);
    case 'fastcore'
        runThresholdsFastCore(figName, bbLetter, modelName, cellLine);
end

end