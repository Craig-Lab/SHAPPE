%% To remove at the end
clearvars
IMCInfo = readtable("SOCPatientsInfo.csv");
inhomogeneityIMC = load('inhomogeneityFinalResults_50.mat');
inhomogeneityIMC = inhomogeneityIMC.inhomogeneityFinalResults;

folder_attractiveRegions_baseline = fullfile(pwd,"AttractiveRegions_DomainEnlargement");
if(~isfolder(folder_attractiveRegions_baseline))
    mkdir(folder_attractiveRegions_baseline)
end

folder_repulsiveRegions_baseline = fullfile(pwd,"RepulsiveRegions_DomainEnlargement");
if(~isfolder(folder_repulsiveRegions_baseline))
    mkdir(folder_repulsiveRegions_baseline)
end

folder_SHAPPEConfig_baseline = fullfile(pwd,"SHAPPEConfigs_DomainEnlargement");
if(~isfolder(folder_SHAPPEConfig_baseline))
    mkdir(folder_SHAPPEConfig_baseline)
end

folder_RandomBaselineConfig_baseline = fullfile(pwd,"RandomBaselineConfigs_DomainEnlargement");
if(~isfolder(folder_RandomBaselineConfig_baseline))
    mkdir(folder_RandomBaselineConfig_baseline)
end

% Loop through sample
elapsedTimes = zeros(1,size(IMCInfo,1));
parfor(np=1:size(IMCInfo,1),20)
    % PP_path = "/Users/morgancraig/Documents/Documents - MacBook Air de wksadm/GBMs/"+...
    %     "ProjectMechanisticReproductionPointPatterns/cellPlacementCode/patientsData/OriginalCSV/"+string(IMCInfo.csv_ID(np))+".csv";
    PP_path = "/NOBACKUP/mongeonb/PointPatternReconstruction_CellPlacement/cellPlacement/patientsData/OriginalCSV/"+string(IMCInfo.csv_ID(np))+".csv";
    Xobs = readmatrix(PP_path);
    Nobs = size(Xobs,1);
    
    %% SECTION 2: (Hyper)Parametrisation
    
    %%% Here, we place cancer cells (type 1) first, then stromal cells (type 3)
    %%% and then immune cells in descending order
    types = 1:6;
    types = types(:).'; % convert to row vector
    G = numel(types);
    NgObs = zeros(1,G);
    for g=1:G
        NgObs(1,g) = sum(Xobs(:,end)==types(g));
    end
    type_proportions = NgObs./Nobs;
    typeOrder = zeros(1,G);
    typeOrder(1) = 1;
    typeOrder(2) = 3;
    typesToRank = setdiff(1:G,[1,3]);
    [~,idx] = sort(type_proportions(typesToRank),'descend');
    for g=(G-numel(typesToRank)+1):G
        typeOrder(g) = typesToRank(idx(g-2));
    end
    
    SHAPPEoptions = struct();
    SHAPPEoptions.typeOrdering = typeOrder;
    
    % Information about reference point pattern (Xobs)
    SHAPPEoptions.RObs = IMCInfo.Radius(np);
    SHAPPEoptions.AreaObs = (2.*(IMCInfo.Radius(np))).^2;
    SHAPPEoptions.ShapeObs = "square";
    SHAPPEoptions.types = types;
    SHAPPEoptions.types_name = ["GBM","CD8","Stroma","M1","M2","M0"];
    
    % Information for inhomogeneity identification test
    SHAPPEoptions.alphaFtest = 0.05; 
    SHAPPEoptions.NRFTest = 1000;
    SHAPPEoptions.NTestLocationsFTest = 10000;
    SHAPPEoptions.A0 = pi.*(50.^2);
    SHAPPEoptions.LowDensityThreshold = 0;
    SHAPPEoptions.bL = 40;
    
    % Information abour random labelling test
    SHAPPEoptions.KFunctionDistances = 15:5:60;
    SHAPPEoptions.BRandomLabelling = 1000;
    SHAPPEoptions.distanceSpanRelationshipSign = [40,50];
    SHAPPEoptions.alphaRandomLabelling = 0.05;
    
    % Information about reconstructions
    SHAPPEoptions.N = 10000;
    SHAPPEoptions.Nreps = 10;
    SHAPPEoptions.shapeSynthetic = "circle";
    SHAPPEoptions.parallelLoops = 0;
    
    folder_attractiveRegions = fullfile(pwd,"AttractiveRegions_DomainEnlargement",string(IMCInfo.csv_ID(np)));
    folder_repulsiveRegions = fullfile(pwd,"RepulsiveRegions_DomainEnlargement",string(IMCInfo.csv_ID(np)));
    folder_SHAPPEConfig = fullfile(pwd,"SHAPPEConfigs_DomainEnlargement",string(IMCInfo.csv_ID(np)));
    folder_RandomBaselineConfig = fullfile(pwd,"RandomBaselineConfigs_DomainEnlargement",string(IMCInfo.csv_ID(np)));
    SHAPPEoptions.outDirectoryRandomBaseline = folder_RandomBaselineConfig;
    SHAPPEoptions.outDirectorySHAPPE = folder_SHAPPEConfig;
    SHAPPEoptions.outDirectoryRepulsiveRegions = folder_repulsiveRegions;
    SHAPPEoptions.outDirectoryAttractiveRegions = folder_attractiveRegions;
    
    % Parameters
    SHAPPEoptions.exclusionDistance = 0.8202;
    SHAPPEoptions.Nmin = 15;
    SHAPPEoptions.ghostDistance = 50;
    SHAPPEoptions.nEnv = 50;
    SHAPPEoptions.bLEnv = 80;
    %%% You can decide to add repulsion from repulsive regions or attraction to
    %%% attractive regions in Step 4 (useRegions = true), or not (useRegions = false)
    SHAPPEoptions.useRegions = true;
    if(np==2)
        SHAPPEoptions.attract = true; %%% Set to false to use repulsive regions and to true to use attractive regions
    else
        SHAPPEoptions.attract = false; %%% Set to false to use repulsive regions and to true to use attractive regions
    end
    SHAPPEoptions.SHAPPEalpha = 0.05;
    
    % Hyperparameters
    SHAPPEoptions.L = 10;
    SHAPPEoptions.a = 0;
    SHAPPEoptions.rW = 1;
    SHAPPEoptions.lambda = 0.01/0.133;
    SHAPPEoptions.phi = 10*10.75;
    SHAPPEoptions.phiMinimum = 10.75;
    SHAPPEoptions.C = "Gaussian";
    SHAPPEoptions.wNoCorr = 0;
    
    % Convergence parameters
    SHAPPEoptions.dMed = 0.1075;
    SHAPPEoptions.d90 = 10.75;
    SHAPPEoptions.dNet = 10.75;
    SHAPPEoptions.nEq = 210;
    SHAPPEoptions.MaxN = 10000;
    
    args = namedargs2cell(SHAPPEoptions);
    
    [SHAPPEoutputs] = SHAPPE(Xobs,args{:});
    elapsedTimes(1,np) = SHAPPEoutputs.elapsedTimeSHAPPE;
end

save("elapsedTimes","elapsedTimes")
