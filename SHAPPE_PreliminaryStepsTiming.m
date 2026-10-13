%% To remove at the end
clearvars
IMCInfo = readtable("SOCPatientsInfo.csv");
inhomogeneityIMC = load('inhomogeneityFinalResults_50.mat');
inhomogeneityIMC = inhomogeneityIMC.inhomogeneityFinalResults;
inhomogeneous = find(inhomogeneityIMC(1,:) & inhomogeneityIMC(2,:));

folder_PreviousCorrelationResults = "/Users/morgancraig/Documents/Documents - MacBook Air de wksadm/GBMs/"+...
    "ProjectMechanisticReproductionPointPatterns/cellPlacementCode/patientsData/CorrelationsAnalysisSquares";

%% Create folders
folder_attractiveRegions = fullfile(pwd,"AttractiveRegions");
folder_repulsiveRegions = fullfile(pwd,"RepulsiveRegions");
folder_SHAPPEConfig = fullfile(pwd,"SHAPPEConfigs");
folder_RandomBaselineConfig = fullfile(pwd,"RandomBaselineConfigs");
if(~isfolder(folder_attractiveRegions))
    mkdir(folder_attractiveRegions)
end

if(~isfolder(folder_repulsiveRegions))
    mkdir(folder_repulsiveRegions)
end

if(~isfolder(folder_SHAPPEConfig))
    mkdir(folder_SHAPPEConfig)
end

if(~isfolder(folder_RandomBaselineConfig))
    mkdir(folder_RandomBaselineConfig)
end

%% Point pattern information and SHAPPE parametrisation
Npatients = size(IMCInfo,1);
inhomogeneityIdentificationTimings = nan(1,Npatients);
wrongInhomogeneityIdentification = [];
randomLabellingTimings = nan(1,Npatients);
wrongRandomLabelling = [];
%%
load("runaug28th.mat")
for np=47:Npatients
    SHAPPEoptions = struct();
    %%% Load your point pattern
    PP_path = "/Users/morgancraig/Documents/Documents - MacBook Air de wksadm/GBMs/"+...
        "ProjectMechanisticReproductionPointPatterns/cellPlacementCode/patientsData/OriginalCSV/"+string(IMCInfo.csv_ID(np))+".csv";
    Xobs = readmatrix(PP_path);
    Nobs = size(Xobs,1);
    
    % Type ordering in point placement: user can define its own ordering
    % preference by passing typeOrder as an array, or by setting it to
    % 'descend' or 'ascend'
    %%% Here, we place cancer cells (type 1) first, then stroma cells (type 3)
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
    SHAPPEoptions.Nreps = 2;
    SHAPPEoptions.shapeSynthetic = "circle";
    SHAPPEoptions.parallelLoops = 1;
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
    SHAPPEoptions.attract = false; %%% Set to false to use repulsive regions and to true to use attractive regions
    SHAPPEoptions.SHAPPEalpha = 0.05;
    
    % Hyperparameters
    SHAPPEoptions.L = 1;
    SHAPPEoptions.a = 0.5;
    SHAPPEoptions.rW = 10;
    SHAPPEoptions.lambda = 0.01/0.133;
    SHAPPEoptions.phi = 2*10.75;
    SHAPPEoptions.phiBaseline = 10.75;
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

    % Check inhomogeneity identification step
    inhomogeneityIdentificationTimings(1,np) = SHAPPEoutputs.HomogeneityIdentificationElapsedTime;
    homogeneousStatus = SHAPPEoutputs.HomogeneityStatus;
    if((homogeneousStatus && ismember(np,inhomogeneous)) || (~homogeneousStatus && ~ismember(np,inhomogeneous)))
        wrongInhomogeneityIdentification = [wrongInhomogeneityIdentification, np];
    end

    % Check random labelling step
    randomLabellingTimings(1,np) = SHAPPEoutputs.RandomLabellingTestElapsedTime;
    % aboveMatrix_optimised = SHAPPEoutputs.aboveMatrix;
    % belowMatrix_optimised = SHAPPEoutputs.belowMatrix;
    % %%% Load old matrices
    % aboveMatrix = loadVariable(fullfile(folder_PreviousCorrelationResults,"above",string(IMCInfo.csv_ID(np))+".mat"));
    % belowMatrix = loadVariable(fullfile(folder_PreviousCorrelationResults,"below",string(IMCInfo.csv_ID(np))+".mat"));
    temp = SHAPPEoutputs.corrMatrix;
    corrMatrix_optimised = zeros(size(temp));
    corrMatrix_optimised(temp<0.05 & ~isnan(temp)) = 1;
    temp = loadVariable(fullfile(folder_PreviousCorrelationResults,string(IMCInfo.csv_ID(np))+".mat"));
    corrMatrix = zeros(size(temp));
    corrMatrix(temp<0.05 & ~isnan(temp)) = 1;

    if(sum(abs(corrMatrix_optimised-corrMatrix),"all")>0)
        wrongRandomLabelling = [wrongRandomLabelling,np];
    end

end

%%
load("runaug28th.mat")

%%
np=5;
SHAPPEoptions = struct();
%%% Load your point pattern
PP_path = "/Users/morgancraig/Documents/Documents - MacBook Air de wksadm/GBMs/"+...
    "ProjectMechanisticReproductionPointPatterns/cellPlacementCode/patientsData/OriginalCSV/"+string(IMCInfo.csv_ID(np))+".csv";
Xobs = readmatrix(PP_path);
Nobs = size(Xobs,1);

% Type ordering in point placement: user can define its own ordering
% preference by passing typeOrder as an array, or by setting it to
% 'descend' or 'ascend'
%%% Here, we place cancer cells (type 1) first, then stroma cells (type 3)
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
SHAPPEoptions.Nreps = 2;
SHAPPEoptions.shapeSynthetic = "circle";
SHAPPEoptions.parallelLoops = 1;
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
SHAPPEoptions.attract = false; %%% Set to false to use repulsive regions and to true to use attractive regions
SHAPPEoptions.SHAPPEalpha = 0.05;

% Hyperparameters
SHAPPEoptions.L = 1;
SHAPPEoptions.a = 0.5;
SHAPPEoptions.rW = 10;
SHAPPEoptions.lambda = 0.01/0.133;
SHAPPEoptions.phi = 2*10.75;
SHAPPEoptions.phiBaseline = 10.75;
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

homogeneousStatus = SHAPPEoutputs.HomogeneityStatus;
% if((homogeneousStatus && ismember(np,inhomogeneous)) || (~homogeneousStatus && ~ismember(np,inhomogeneous)))
%     wrongInhomogeneityIdentification = [wrongInhomogeneityIdentification, np];
% end

% Check random labelling step
temp = SHAPPEoutputs.corrMatrix;
corrMatrix_optimised = zeros(size(temp));
corrMatrix_optimised(temp<0.05 & ~isnan(temp)) = 1;
temp = loadVariable(fullfile(folder_PreviousCorrelationResults,string(IMCInfo.csv_ID(np))+".mat"));
corrMatrix = zeros(size(temp));
corrMatrix(temp<0.05 & ~isnan(temp)) = 1;

corrMatrix = loadVariable(fullfile(folder_PreviousCorrelationResults,string(IMCInfo.csv_ID(np))+".mat"));
aboveMatrix_optimised = SHAPPEoutputs.aboveMatrix;
belowMatrix_optimised = SHAPPEoutputs.belowMatrix;
%%% Load old matrices
aboveMatrix = loadVariable(fullfile(folder_PreviousCorrelationResults,"above",string(IMCInfo.csv_ID(np))+".mat"));
belowMatrix = loadVariable(fullfile(folder_PreviousCorrelationResults,"below",string(IMCInfo.csv_ID(np))+".mat"));
