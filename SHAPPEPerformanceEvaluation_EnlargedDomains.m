close all
clearvars



%% Retrieve paths to folder
folder_code = pwd;
folder_SHAPPEconfig = fullfile(folder_code,"SHAPPEConfigs_DomainEnlargement");
folder_patientsData = "/NOBACKUP/mongeonb/PointPatternReconstruction_CellPlacement/cellPlacement/patientsData";
folder_originalCSV = "/NOBACKUP/mongeonb/PointPatternReconstruction_CellPlacement/cellPlacement/patientsData/OriginalCSV/";
folder_randomBaseline = fullfile(folder_code,"RandomBaselineConfigs_DomainEnlargement");

%% Constants
%Indexes for cell types in original csv file
csv.CD8 = 2;
csv.GBM = 1;
csv.Stroma = 3;
csv.M0 = 6;
csv.M1 = 4;
csv.M2 = 5;
csv_cellTypes =[csv.GBM, csv.CD8, csv.Stroma, csv.M0, csv.M1, csv.M2];
celltypes_name = ["GBM","CD8","Stroma","M0","M1","M2"];
Ncelltypes = length(celltypes_name);

%Number of cells in augmented configurations
N=10000;

%Cell radius
radius.GBM = 10.75;
radius.CD8 = 3.6;
radius.Stroma = 7.48;
radius.M1 = 5;
radius.M2 = 5;
radius.M0 = 5;
radii = zeros(1,Ncelltypes);
for t=1:Ncelltypes
    ct = celltypes_name(t);
    radii(csv.(ct)) = radius.(ct);
end % End of for t=1:Ncelltypes

%% Load file with patients' information
patientInfo_together = readtable(fullfile(folder_code,'SOCPatientsInfo.csv'));
Npatients = size(patientInfo_together,1);

%% SECTION 1: Calculate metrics

%%%%% Set constants
coldThreshold = 0;
cancerCellType = csv.GBM;
rneigh = 50;
rL = [20 40 60];
minCellCount = 15;
nReplicates = 10;
bLanalysis = 40;
sigma_micron = 10./2.355;
h_kernel = 100;

fileLBs = fullfile(folder_patientsData,"RipleyLBs_markBootstrap.mat");
load(fileLBs)
fileUBs = fullfile(folder_patientsData,"RipleyUBs_markBootstrap.mat");
load(fileUBs)
fileFonctionK = fullfile(folder_patientsData,"RipleyK_Patients.mat");
load(fileFonctionK)
load(fullfile(folder_patientsData,"inhomogeneityFinalResults_50.mat"))

%%%%% Set jobs
patients = (1:Npatients)';
allJobs = [patients,ones(size(patients));...
    patients,2.*ones(size(patients))];
nJobs = size(allJobs,1);

%%%%% Define metrics
metricsToCalculate = struct();
metricsToCalculate.Ldiff = true;
metricsToCalculate.OV_NN = true;
metricsToCalculate.ExposureMatrices = true;
metricsToCalculate.DensityMapDifference = true;
metricsToCalculate.coldRegionsMetrics = false;

resultsRipleyL1 = nan(nJobs,1);
resultsRipleyL2 = nan(nJobs,1);
resultsRipleyL3 = nan(nJobs,1);
% Out of CIs
resultsRipleyL1_exceedance = nan(nJobs,1);
resultsRipleyL2_exceedance = nan(nJobs,1);
resultsRipleyL3_exceedance = nan(nJobs,1);

ovNNDiff = nan(nJobs,1);
DMDiff = nan(nJobs,1);
resultsscoreMu = nan(nJobs,1);
resultsscorePi = nan(nJobs,1);
resultsscoreE = nan(nJobs,1);

done = false(nJobs,1); 

cd(folder_code)

fileResults = fullfile(folder_code,'SHAPPEPerformance_EnlargedDomains.mat');

if(isfile(fileResults))
    load(fileResults,"resultsRipleyL1","resultsRipleyL2","resultsRipleyL3",...
        "resultsRipleyL1_exceedance","resultsRipleyL2_exceedance","resultsRipleyL3_exceedance",...
        "ovNNDiff","DMDiff","resultsscoreMu","resultsscorePi","resultsscoreE","done")
end

% How many jobs left to do
uncompletedJobs = find(~done);
NuncompledJobs = length(uncompletedJobs);
uncompletedJobsInfo = allJobs(uncompletedJobs,:);

% Loop through jobs
M = 25;
nJobsPerPartition = NuncompledJobs;
nPartition = ceil(NuncompledJobs/nJobsPerPartition);

for p=1:nPartition
    start_index_loop = (p-1).*nJobsPerPartition+1;
    end_index_loop = min(p.*nJobsPerPartition,NuncompledJobs);
    disp("Starting partition #"+string(p)+".")

    jobsInPartition = uncompletedJobs(start_index_loop:end_index_loop);
    jobsInPartition = jobsInPartition(:);

    temp_resultsRipleyL1 = nan(length(jobsInPartition),1);
    temp_resultsRipleyL2 = nan(length(jobsInPartition),1);
    temp_resultsRipleyL3 = nan(length(jobsInPartition),1);
    temp_resultsRipleyL1_exceedance = nan(length(jobsInPartition),1);
    temp_resultsRipleyL2_exceedance = nan(length(jobsInPartition),1);
    temp_resultsRipleyL3_exceedance = nan(length(jobsInPartition),1);
    temp_ovNNDiff = nan(length(jobsInPartition),1);
    temp_DMDiff = nan(length(jobsInPartition),1);
    temp_resultsscoreMu = nan(length(jobsInPartition),1);
    temp_resultsscorePi = nan(length(jobsInPartition),1);
    temp_resultsscoreE = nan(length(jobsInPartition),1);
    temp_done = false(length(jobsInPartition),1);

    parfor (jj = 1:length(jobsInPartition), M)
        jobIndex = jobsInPartition(jj);
        currentJob = allJobs(jobIndex,:);
        np = currentJob(1,1);
        scenario = currentJob(1,2);

        current_px = patientInfo_together(np,:);
        patientName = string(current_px.csv_ID);

        px_csv = strcat(patientName,".csv");
        filenameIMC = fullfile(folder_originalCSV,px_csv);
        originalIMCConfiguration = readmatrix(filenameIMC);
        radiusIMC = current_px.Radius;
        radiusABM = current_px.synthetic_radius_circle;
        cellProportions = [current_px.proportion_GBM current_px.proportion_Stroma ...
            current_px.proportion_CD8 current_px.proportion_M1 ...
            current_px.proportion_M2 current_px.proportion_M0];
        NcellsIMC = zeros(1,6);
        for nt=1:length(csv_cellTypes)
            ct = csv_cellTypes(nt);
            NcellsIMC(ct) = sum(originalIMCConfiguration(:,end)==ct);
        end

        RipleyExcess_CI = zeros(1,length(rL));
        RipleyDifference_raw = zeros(1,length(rL));
        ripleyLweights = ones(size(cellProportions));
        LB_CurrentPatient= reshape(LBRipley_allPatients(np,:,:,:),[Ncelltypes,Ncelltypes,length(rL)]);
        UB_CurrentPatient = reshape(UBRipley_allPatients(np,:,:,:),[Ncelltypes,Ncelltypes,length(rL)]);
        OV_NN = 0;
        scoreMu = 0;
        scorePi = 0;
        scoreE = 0;
        densityMapDifference = 0;

        % Initialise structure with IMC sample information
        informationIMC = struct();
        informationIMC.configuration = originalIMCConfiguration;
        informationIMC.Ncells = NcellsIMC;
        informationIMC.radius = radiusIMC;
        informationIMC.binLength = bLanalysis;
        informationIMC.sigma_micron = sigma_micron;
        informationIMC.DensityThreshold = coldThreshold;
        informationIMC.InhomogeneityStatus = inhomogeneityFinalResults(1,np)>0 && inhomogeneityFinalResults(2,np)>0;
        informationIMC.RipleyK = reshape(Rippley_allPatient(np,:,:,:),[Ncelltypes,Ncelltypes,length(rL)]);
        informationIMC.h_kernel = h_kernel;
        informationIMC.shape = "square";

        if(scenario==1)
            folder_replicates_px = fullfile(folder_SHAPPEconfig,patientName);
            nCompletedReplicates = nReplicates;
        elseif(scenario==2)
            folder_replicates_px = fullfile(folder_randomBaseline,patientName);
            nCompletedReplicates = nReplicates;
        end % End of if(scenario==1) ... elseif(scenario==2) end

        % Store the number of valid replicates in case something goes wrong
        num_replicates = 0;
        NcellsABM = zeros(1,6);
        if(nCompletedReplicates == nReplicates)
            
            for nr = 1:nReplicates
                filenameConfiguration = fullfile(folder_replicates_px,...
                    strcat("rep",num2str(nr),".csv"));
                fileInfoConfiguration = dir(filenameConfiguration);
                if(isfile(filenameConfiguration) && ~isempty(fileInfoConfiguration) && fileInfoConfiguration.bytes > 0)
                    num_replicates = num_replicates+1;
                    configuration = readmatrix(filenameConfiguration);
                    NcellsABM = NcellsIMC;
                    for nt=1:length(csv_cellTypes)
                        ct = csv_cellTypes(nt);
                        NcellsABM(ct) = sum(configuration(:,end)==ct);
                    end
                    informationConfiguration = struct();
                    informationConfiguration.configuration = configuration;
                    informationConfiguration.Ncells = NcellsABM;
                    informationConfiguration.radius = radiusABM;
                    informationConfiguration.shape = "circle";

                    metricsConfiguration = SHAPPE_Metrics(informationIMC,informationConfiguration,rL,ripleyLweights,...
                        LB_CurrentPatient,UB_CurrentPatient,minCellCount,csv_cellTypes,rneigh,metricsToCalculate);

                    % Ripley's L raw difference and exceedance out of CIs (0 if inside CI)
                    RipleyDifference_raw= RipleyDifference_raw+metricsConfiguration.RipleyDifference_raw;
                    RipleyExcess_CI = RipleyExcess_CI+metricsConfiguration.RipleyExcess_CI;
                    % Overlap with nearest neighbor distribution
                    OV_NN = OV_NN+metricsConfiguration.OV_NN;
                    %%%%% Density Map Difference
                    densityMapDifference = densityMapDifference+metricsConfiguration.densityMapDifference;
                    %%%%% Exposure statistics
                    scoreMu = scoreMu+metricsConfiguration.score_mu;
                    scorePi = scorePi+metricsConfiguration.score_pi;
                    scoreE = scoreE+metricsConfiguration.score_e;
                end % End of if(isfile(filenameConfiguration) && ~isempty(fileInfoConfiguration) && fileInfoConfiguration.bytes > 0)
            end % End of for nr = 1:nReplicates
            if(num_replicates>0)
                temp_done(jj) = true;
                % Average over num_replicates
                temp_resultsRipleyL1(jj) = RipleyDifference_raw(1)/num_replicates;
                temp_resultsRipleyL2(jj) = RipleyDifference_raw(2)/num_replicates;
                temp_resultsRipleyL3(jj) = RipleyDifference_raw(3)/num_replicates;
                temp_resultsRipleyL1_exceedance(jj) = RipleyExcess_CI(1)/num_replicates;
                temp_resultsRipleyL2_exceedance(jj) = RipleyExcess_CI(2)/num_replicates;
                temp_resultsRipleyL3_exceedance(jj) = RipleyExcess_CI(3)/num_replicates;
                temp_resultsscoreMu(jj) = scoreMu/num_replicates;
                temp_resultsscorePi(jj) = scorePi/num_replicates;
                temp_resultsscoreE(jj) = scoreE/num_replicates;
                temp_ovNNDiff(jj) = OV_NN/num_replicates;
                temp_DMDiff(jj) = densityMapDifference/num_replicates;
            end % End of if(num_replicates>0)
        end % End of if(nCompletedReplicates == nReplicates)
    end % End of parfor(jj=start_index_loop:end_index_loop,M)
    resultsRipleyL1(jobsInPartition) = temp_resultsRipleyL1;
    resultsRipleyL2(jobsInPartition) = temp_resultsRipleyL2;
    resultsRipleyL3(jobsInPartition) = temp_resultsRipleyL3;
    resultsRipleyL1_exceedance(jobsInPartition) = temp_resultsRipleyL1_exceedance;
    resultsRipleyL2_exceedance(jobsInPartition) = temp_resultsRipleyL2_exceedance;
    resultsRipleyL3_exceedance(jobsInPartition) = temp_resultsRipleyL3_exceedance;
    resultsscoreMu(jobsInPartition) = temp_resultsscoreMu;
    resultsscorePi(jobsInPartition) = temp_resultsscorePi;
    resultsscoreE(jobsInPartition) = temp_resultsscoreE;
    ovNNDiff(jobsInPartition) = temp_ovNNDiff;
    DMDiff(jobsInPartition) = temp_DMDiff;
    done(jobsInPartition) = temp_done;

    save(fileResults)
end




