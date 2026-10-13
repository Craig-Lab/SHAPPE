close all
clearvars
addpath(fileparts(pwd))
%% Load results

load('SHAPPE_HyperparameterOptimisation_SameSize.mat')
% load('SHAPPE_HyperparameterOptimisation_DomainIncrease.mat')
% load('SHAPPEPerformance_EnlargedDomains.mat')
%%
inhomogeneityIdentification = loadVariable("inhomogeneityFinalResults_50.mat");
InhomogeneousPatients = find(inhomogeneityIdentification(1,:)>0 & inhomogeneityIdentification(2,:)>0);
HomogeneousPatients = find(inhomogeneityIdentification(1,:)==0 | inhomogeneityIdentification(2,:)==0);
set_colors;

patientInfo = readtable("SOCPatientsInfo.csv");
Npatients = size(patientInfo,1);

cancerCellType = 1;
celltypes_name = ["GBM","CD8","Stroma","M0","M1","M2"];
csv.GBM = 1;
csv.CD8 = 2;
csv.Stroma = 3;
csv.M1 = 4;
csv.M2 = 5;
csv.M0 = 6;
csv_cellTypes = [csv.GBM csv.CD8 csv.Stroma csv.M0 csv.M1 csv.M2];
Ncelltypes = numel(csv_cellTypes);

radius.GBM = 10.75;
radius.CD8 = 3.6;
radius.Stroma = 7.48;
radius.M1 = 5;
radius.M2 = 5;
radius.M0 = 5;
radii = [radius.GBM radius.CD8 radius.Stroma radius.M0 radius.M1 radius.M2];

col_homogeneous = oranges(5,:);
col_inhomogeneous = purples(2,:);

Nreplicates = 10;

%% set folders
folders.originalCSV = "OriginalCSV";
folders.ECSRConfigurations = "ECSR_SameSize_Square";
folders.SHAPPEOptimisedConfig = "SHAPPEConfig_SameSize_Square";
folders.SHAPPESimplifiedConfig = "SHAPPESimplifiedConfig_SameSize_Square";

folders.figures = "/Users/morgancraig/Documents/Documents - MacBook Air de wksadm/"+...
    "GBMs/ProjectMechanisticReproductionPointPatterns/figures";
folders.MAINTEXTFigures = fullfile(folders.figures,"MAINTEXT");
folders.ABMFigures = fullfile(folders.MAINTEXTFigures,"ABM_NewABM");
folders.DensityFigures = fullfile(folders.MAINTEXTFigures,"DensityPatternReproduction");
folders.ClusteringFigures = fullfile(folders.MAINTEXTFigures,"FIGURE_Clustering_RandomVsHooke");
% folders.FrameworkFigures = fullfile(folders.MAINTEXTFigures,"Framework");
% folders.InhomogeneityIdentificationProcessFigures = fullfile(folders.MAINTEXTFigures,"InhomogeneityIdentificationProcess");
folders.SIFigures = fullfile(folders.figures,"SI_Figures");
folders.ExposureFigures = fullfile(folders.SIFigures,"Exposure");
folders.MetricsFigures = fullfile(folders.SIFigures,"Metrics");
folders.HyperparameterOptimisationFigures = fullfile(folders.SIFigures,"HyperparameterOptimisation");
folders.InsideBoundsFigures = fullfile(folders.SIFigures,"InsideBounds");
folders.ExecutionConvergenceFigures = fullfile(folders.SIFigures,"ExecutionConvergence");
% folders.BoundaryConditionsFigures = fullfile(folders.SIFigures,"BoundaryConditions");
% folders.CellCyclesABMFigures = fullfile(folders.SIFigures,"CellCyclesABM");
% folders.CSRFigures = fullfile(folders.SIFigures,"CSR");
% folders.ECSRFigures = fullfile(folders.SIFigures,"ECSR");
% folders.GhostPointsFigures = fullfile(folders.SIFigures,"GhostPoints");
% folders.InhomogeneityIdentificationProcess_ResultsFigures = fullfile(folders.SIFigures,"InhomogeneityIdentificationProcess_Results");
% folders.InverseTransformMethodFigures = fullfile(folders.SIFigures,"InverseTransformMethod");
% folders.BoundaryConditionsFigures = fullfile(folders.SIFigures,"BoundaryConditions");
% folders.RandomLabellingFigures = fullfile(folders.SIFigures,"RandomLabelling");
% folders.voidAcceptanceProbabilityFigures = fullfile(folders.SIFigures,"voidAcceptanceProbability");

fields = string(fieldnames(folders));
for ff=1:numel(fields)
    f = folders.(fields(ff));
    if(~isfolder(f))
        mkdir(f)
    end
end

%% Reshape results
SHAPPEIndex = 1:1:Npatients.*nComb;
%%% Homogeneous Random: All CSR
HomogeneousRandomIndex = (SHAPPEIndex(end)+1):1:(SHAPPEIndex(end)+Npatients);
%%% Heterogeneous Random: CSR for homogeneous samples and ECSR for
%%% inhomogeneous samples
HeterogeneousRandomIndex = (HomogeneousRandomIndex(end)+1):1:(HomogeneousRandomIndex(end)+Npatients);

%%%%% SHAPPE
%%% Ripley's H function - Excess out of CI
% SHAPPE
metricsSameSize.H20_CIexcess.SHAPPE = reshape(RipleyH_20_OutOfCI(SHAPPEIndex),[nComb Npatients]);
metricsSameSize.H40_CIexcess.SHAPPE = reshape(RipleyH_40_OutOfCI(SHAPPEIndex),[nComb Npatients]);
metricsSameSize.H60_CIexcess.SHAPPE = reshape(RipleyH_60_OutOfCI(SHAPPEIndex),[nComb Npatients]);
% CSR
metricsSameSize.H20_CIexcess.CSR = reshape(RipleyH_20_OutOfCI(HomogeneousRandomIndex),[1 Npatients]);
metricsSameSize.H40_CIexcess.CSR = reshape(RipleyH_40_OutOfCI(HomogeneousRandomIndex),[1 Npatients]);
metricsSameSize.H60_CIexcess.CSR = reshape(RipleyH_60_OutOfCI(HomogeneousRandomIndex),[1 Npatients]);
% ECSR/CSR
metricsSameSize.H20_CIexcess.ECSR = reshape(RipleyH_20_OutOfCI(HeterogeneousRandomIndex),[1 Npatients]);
metricsSameSize.H40_CIexcess.ECSR = reshape(RipleyH_40_OutOfCI(HeterogeneousRandomIndex),[1 Npatients]);
metricsSameSize.H60_CIexcess.ECSR = reshape(RipleyH_60_OutOfCI(HeterogeneousRandomIndex),[1 Npatients]);
%%% Ripley's H function - Raw Differences
% SHAPPE
metricsSameSize.H20Diff.SHAPPE = reshape(RipleyH_20_RawDifference(SHAPPEIndex),[nComb Npatients]);
metricsSameSize.H40Diff.SHAPPE = reshape(RipleyH_40_RawDifference(SHAPPEIndex),[nComb Npatients]);
metricsSameSize.H60Diff.SHAPPE = reshape(RipleyH_60_RawDifference(SHAPPEIndex),[nComb Npatients]);
% CSR
metricsSameSize.H20Diff.CSR = reshape(RipleyH_20_RawDifference(HomogeneousRandomIndex),[1 Npatients]);
metricsSameSize.H40Diff.CSR = reshape(RipleyH_40_RawDifference(HomogeneousRandomIndex),[1 Npatients]);
metricsSameSize.H60Diff.CSR = reshape(RipleyH_60_RawDifference(HomogeneousRandomIndex),[1 Npatients]);
% ECSR/CSR
metricsSameSize.H20Diff.ECSR = reshape(RipleyH_20_RawDifference(HeterogeneousRandomIndex),[1 Npatients]);
metricsSameSize.H40Diff.ECSR = reshape(RipleyH_40_RawDifference(HeterogeneousRandomIndex),[1 Npatients]);
metricsSameSize.H60Diff.ECSR = reshape(RipleyH_60_RawDifference(HeterogeneousRandomIndex),[1 Npatients]);
%%% Density Map Differences
% SHAPPE
metricsSameSize.DMDiff.SHAPPE = reshape(DMDiff(SHAPPEIndex,1),[nComb Npatients]);
% CSR
metricsSameSize.DMDiff.CSR = reshape(DMDiff(HomogeneousRandomIndex,1),[1 Npatients]);
% ECSR/CSR
metricsSameSize.DMDiff.ECSR = reshape(DMDiff(HeterogeneousRandomIndex,1),[1 Npatients]);
%%% NN Difference
% SHAPPE
metricsSameSize.NNDiff.SHAPPE = reshape(NNDiff(SHAPPEIndex,1),[nComb Npatients]);
% CSR
metricsSameSize.NNDiff.CSR = reshape(NNDiff(HomogeneousRandomIndex,1),[1 Npatients]);
% ECSR/CSR
metricsSameSize.NNDiff.ECSR = reshape(NNDiff(HeterogeneousRandomIndex,1),[1 Npatients]);
%%% Exposure Status
% SHAPPE
metricsSameSize.Sigma.SHAPPE = reshape(ExposureStatus(SHAPPEIndex,1),[nComb Npatients]);
% CSR
metricsSameSize.Sigma.CSR = reshape(ExposureStatus(HomogeneousRandomIndex,1),[1 Npatients]);
% ECSR/CSR
metricsSameSize.Sigma.ECSR = reshape(ExposureStatus(HeterogeneousRandomIndex,1),[1 Npatients]);
%%% Exposure Levels
% SHAPPE
metricsSameSize.Pi.SHAPPE = reshape(ExposureLevels(SHAPPEIndex,1),[nComb Npatients]);
% CSR
metricsSameSize.Pi.CSR = reshape(ExposureLevels(HomogeneousRandomIndex,1),[1 Npatients]);
% ECSR/CSR
metricsSameSize.Pi.ECSR = reshape(ExposureLevels(HeterogeneousRandomIndex,1),[1 Npatients]);

%% Define Objective Functions

%%% Objective function 2: H(20) (Raw difference)
F2 = metricsSameSize.H20Diff.SHAPPE;
F2_ub = max(F2(:));
F2_lb = 0;
F2_normd = (F2-F2_lb)./(F2_ub-F2_lb);
F2_avg = mean(F2_normd,2,"omitnan");

%%% Objective function 3: H(40) (Raw difference)
F3 = metricsSameSize.H40Diff.SHAPPE;
F3_ub = max(F3(:));
F3_lb = 0;
F3_normd = (F3-F3_lb)./(F3_ub-F3_lb);
F3_avg = mean(F3_normd,2,"omitnan");

%%% Objective function 4: H(60) (Raw difference)
F4 = metricsSameSize.H60Diff.SHAPPE;
F4_ub = max(F4(:));
F4_lb = 0;
F4_normd = (F4-F4_lb)./(F4_ub-F4_lb);
F4_avg = mean(F4_normd,2,"omitnan");
%% Pareto Optimisation

minimisationResults = struct();

%%% Ripley's H functions at r=20, r=40 and r=60
%%%%% Homogeneous and Inhomogeneous
scoresdata = [F2_avg,F3_avg,F4_avg];
[ParetoMembership,Pareto_Member_value] = find_pareto_frontier(scoresdata);
idxParetoMembership = find(ParetoMembership);
L2normParetoFrontier = vecnorm(Pareto_Member_value,2,2);
[~,idxMinParetoFrontier] = min(L2normParetoFrontier);
idxOptimised = idxParetoMembership(idxMinParetoFrontier);
minimisationResults.Function2.Both.idxParetoMembership = idxParetoMembership;
minimisationResults.Function2.Both.idxOptimised = idxOptimised;

%% Put results together
idxOptimised = minimisationResults.Function2.Both.idxOptimised;
idxH20_CIexcess = 1;
idxH40_CIexcess = 2;
idxH60_CIexcess = 3;
idxH20_rawDiff = 4;
idxH40_rawDiff = 5;
idxH60_rawDiff = 6;
idxDMDiff = 7;
idxNNDiff = 8;
idxSigma = 9;
idxPi = 10;

%%% For all point patterns
resultsOptimised = [metricsSameSize.H20_CIexcess.SHAPPE(idxOptimised,:);...
    metricsSameSize.H40_CIexcess.SHAPPE(idxOptimised,:);...
    metricsSameSize.H60_CIexcess.SHAPPE(idxOptimised,:);...
    metricsSameSize.H20Diff.SHAPPE(idxOptimised,:);...
    metricsSameSize.H40Diff.SHAPPE(idxOptimised,:);...
    metricsSameSize.H60Diff.SHAPPE(idxOptimised,:);...
    metricsSameSize.DMDiff.SHAPPE(idxOptimised,:);...
    metricsSameSize.NNDiff.SHAPPE(idxOptimised,:);...
    metricsSameSize.Sigma.SHAPPE(idxOptimised,:);...
    metricsSameSize.Pi.SHAPPE(idxOptimised,:)];
resultsCSR = [metricsSameSize.H20_CIexcess.CSR;...
    metricsSameSize.H40_CIexcess.CSR;...
    metricsSameSize.H60_CIexcess.CSR;...
    metricsSameSize.H20Diff.CSR;...
    metricsSameSize.H40Diff.CSR;...
    metricsSameSize.H60Diff.CSR;...
    metricsSameSize.DMDiff.CSR;...
    metricsSameSize.NNDiff.CSR;...
    metricsSameSize.Sigma.CSR;...
    metricsSameSize.Pi.CSR];
resultsECSR = [metricsSameSize.H20_CIexcess.ECSR;...
    metricsSameSize.H40_CIexcess.ECSR;...
    metricsSameSize.H60_CIexcess.ECSR;...
    metricsSameSize.H20Diff.ECSR;...
    metricsSameSize.H40Diff.ECSR;...
    metricsSameSize.H60Diff.ECSR;...
    metricsSameSize.DMDiff.ECSR;...
    metricsSameSize.NNDiff.ECSR;...
    metricsSameSize.Sigma.ECSR;...
    metricsSameSize.Pi.ECSR];

%%% Without additional environmental constraints in Step 4
resultsSimplified = load('SHAPPE_HyperparameterOptimisation_Simplified.mat','resultsOptimised');
resultsSimplified = resultsSimplified.resultsOptimised;

%%% Larger configurations
resultsDomainIncrease = load('SHAPPE_HyperparameterOptimisation_DomainIncrease_new.mat','resultsOptimised');
resultsDomainIncrease = resultsDomainIncrease.resultsOptimised;

%% Sort patients
scoresH = mean([resultsOptimised(idxH20_rawDiff,:);...
    resultsOptimised(idxH40_rawDiff,:);...
    resultsOptimised(idxH60_rawDiff,:)]);
[~,sortIdx] = sort(scoresH,'descend');

%% Point pattern specific optimisation
resultsPPSpecific = zeros(size(resultsOptimised));
optimisedCombNumber_PPSpecific = zeros(1,Npatients);
optimisedSets_PPSpecific = zeros(size(paramComb,2),Npatients);

for np=1:Npatients
    scoresdata = [F2_normd(:,np),F3_normd(:,np),F4_normd(:,np)];
    [ParetoMembership,Pareto_Member_value] = find_pareto_frontier(scoresdata);
    idxParetoMembership = find(ParetoMembership);
    L2normParetoFrontier = vecnorm(Pareto_Member_value,2,2);
    [~,idxMinParetoFrontier] = min(L2normParetoFrontier);
    idxOptimised = idxParetoMembership(idxMinParetoFrontier);
    optimisedCombNumber_PPSpecific(1,np) = idxOptimised;
    optimisedSets_PPSpecific(:,np) = paramComb(idxOptimised,:)';

    resultsPPSpecific(:,np) = [metricsSameSize.H20_CIexcess.SHAPPE(idxOptimised,np);...
        metricsSameSize.H40_CIexcess.SHAPPE(idxOptimised,np);...
        metricsSameSize.H60_CIexcess.SHAPPE(idxOptimised,np);...
        metricsSameSize.H20Diff.SHAPPE(idxOptimised,np);...
        metricsSameSize.H40Diff.SHAPPE(idxOptimised,np);...
        metricsSameSize.H60Diff.SHAPPE(idxOptimised,np);...
        metricsSameSize.DMDiff.SHAPPE(idxOptimised,np);...
        metricsSameSize.NNDiff.SHAPPE(idxOptimised,np);...
        metricsSameSize.Sigma.SHAPPE(idxOptimised,np);...
        metricsSameSize.Pi.SHAPPE(idxOptimised,np)];
end

%% Save results

resultsOptimisation.PPspecific = resultsPPSpecific;
resultsOptimisation.All = resultsOptimised;
resultsOptimisation.Random = resultsECSR;
resultsOptimisation.Additional = resultsSimplified;
resultsOptimisation.DomainIncrease = resultsDomainIncrease;
save("resultsOptimisation.mat",'resultsOptimisation')

%% Create a csv for hyperparameter importance analysis

%%% Columns: Copula - phi - rw - L - wNoCorr - alpha
paramComb_realValues = paramComb;
paramComb_realValues(:,5) = paramComb(:,5)./10;
hyperparametersResultTable = array2table(paramComb_realValues,...
    VariableNames=["Copula","phi" "rw","L","wNoCorr","alpha"]);

%%% Add a column for the objective function
scoresdata = [F2_avg,F3_avg,F4_avg];
combinedScore = vecnorm(scoresdata,2,2);
hyperparametersResultTable_combinedScores = addvars(hyperparametersResultTable,combinedScore,NewVariableNames="RipleyRawDifference_final");

%%% Save results
writetable(hyperparametersResultTable_combinedScores,"hyperparametersResultTable_combinedScores_final.csv")

%%% Do the same with multiple objective functions
hyperparametersResultTable_multipleScores = addvars(hyperparametersResultTable,F2_avg,F3_avg,F4_avg,NewVariableNames=["H20","H40","H60"]);

%%% Save results
writetable(hyperparametersResultTable_multipleScores,"hyperparametersResultTable_multipleScores_final.csv")

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLUSTERING RESULT %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NN-Distribution
exportFigureNN = false;

figNN_BaselineComparison = figure;
hold on
yR = resultsECSR(idxNNDiff,sortIdx);
yH = resultsOptimised(idxNNDiff,sortIdx);
den = yR;
den(den==0) = 1;
y = (yR-yH)./den;
plotIQR(y,figNN_BaselineComparison,gray);
b = bar(y,'FaceColor','flat','BarWidth',1);
for np=1:Npatients
    if(ismember(sortIdx(np),HomogeneousPatients))
        col = col_homogeneous;
    else
        col = col_inhomogeneous;
    end
    b.CData(np,:) = col;
end
yline(0)
xticks([])
xlabel("IMC Point Patterns")
ylabel("\Delta e_{NN}/e_{NN}^{(ECSR)}")
set(gca,'FontSize',axis_fontsize,'FontName',axis_fontname,...
    'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
height = 3.66;
width  = 13.64;
set(figNN_BaselineComparison, 'Units',  'Centimeters')
set(figNN_BaselineComparison, 'Position',  [0, 0,width,height])
set(figNN_BaselineComparison,'Units','Inches');
pos = get(figNN_BaselineComparison,'Position');
set(figNN_BaselineComparison,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
printFig(figNN_BaselineComparison)
if(exportFigureNN)
    exportgraphics(figNN_BaselineComparison,fullfile(folders.ClusteringFigures,"eNN.pdf"),...
        'Resolution',600,'BackgroundColor', 'none', 'ContentType', 'vector',...
        'Width',width, 'Height',height, 'Units','centimeters', 'Padding','figure')
end


%% Ripley H function: Optimisation on all point patterns
%%% Figure 1: Raw H(20), H(40) and H(60)
exportFigureH = false;
scoreIndex = [idxH20_rawDiff,idxH40_rawDiff,idxH60_rawDiff];
scoreNames = ["\Delta e_{H}(20)","\Delta e_{H}(40)","\Delta e_{H}(60)"];
figRipleyH = figure;
t = tiledlayout(length(scoreIndex),1);
allaxes = [];
ylims = zeros(3,2);
for ii=1:length(scoreIndex)
    ax = nexttile(ii);
    allaxes = [allaxes,ax];
    hold on
    yR = resultsECSR(scoreIndex(ii),sortIdx);
    yH = resultsOptimised(scoreIndex(ii),sortIdx);
    den = yR;
    den(den==0) = 1;
    y = (yR-yH);
    plotIQR(y,figRipleyH,gray);
    b = bar(y,'FaceColor','flat','BarWidth',1);
    for np=1:Npatients
        if(ismember(sortIdx(np),HomogeneousPatients))
            col = col_homogeneous;
        else
            col = col_inhomogeneous;
        end
        b.CData(np,:) = col;
    end
    yline(0)
    ylabel(scoreNames(ii))
    xticks([])
    if(ii==length(scoreIndex))
        xlabel("IMC Point Patterns")
    end
    yl = ylim;
    ylims(ii,:) = yl;
    ylim(yl)
end

set(allaxes,'FontSize',axis_fontsize,'FontName',axis_fontname,...
    'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
height = 6.91;
width  = 6.81;
set(figRipleyH, 'Units',  'Centimeters')
set(figRipleyH, 'Position',  [0, 0,width,height])
set(figRipleyH,'Units','Inches');
pos = get(figRipleyH,'Position');
set(figRipleyH,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
printFig(figRipleyH)
if(exportFigureH)
    exportgraphics(figRipleyH, fullfile(folders.ClusteringFigures,"RipleyH.pdf"),...
        'ContentType','vector', 'BackgroundColor','none', ...
        'Width',width, 'Height',height, 'Units','centimeters', 'Padding','figure')
end

%% Ripley H function: Point Pattern Specific Optimisation
%%% Figure 1: Raw H(20), H(40) and H(60)
exportFigureHPPSpecific= false;
scoreIndex = [idxH20_rawDiff,idxH40_rawDiff,idxH60_rawDiff];
scoreNames = ["\Delta e_{H}(20)","\Delta e_{H}(40)","\Delta e_{H}(60)"];
figH_PPSpecific = figure;
t = tiledlayout(length(scoreIndex),1);
allaxes = [];
for ii=1:length(scoreIndex)
    ax = nexttile(ii);
    allaxes = [allaxes,ax];
    hold on
    yR = resultsECSR(scoreIndex(ii),sortIdx);
    yH = resultsPPSpecific(scoreIndex(ii),sortIdx);
    den = yR;
    den(den==0) = 1;
    y = (yR-yH);
    plotIQR(y,figH_PPSpecific,gray);
    b = bar(y,'FaceColor','flat','BarWidth',1);
    for np=1:Npatients
        if(ismember(sortIdx(np),HomogeneousPatients))
            col = col_homogeneous;
        else
            col = col_inhomogeneous;
        end
        b.CData(np,:) = col;
    end
    yline(0)
    ylabel(scoreNames(ii))
    xticks([])
    if(ii==length(scoreIndex))
        xlabel("IMC Point Patterns")
    end
    current_yl = ylim;
    ylim([min(current_yl(1),ylims(ii,1)) max(current_yl(2),ylims(ii,2))])
end
set(allaxes,'FontSize',axis_fontsize,'FontName',axis_fontname,...
    'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
height = 6.91;
width  = 6.81;
set(figH_PPSpecific, 'Units',  'Centimeters')
set(figH_PPSpecific, 'Position',  [0, 0,width,height])
set(figH_PPSpecific,'Units','Inches');
pos = get(figH_PPSpecific,'Position');
set(figH_PPSpecific,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
printFig(figH_PPSpecific)
if(exportFigureHPPSpecific)
    exportgraphics(figH_PPSpecific, fullfile(folders.ClusteringFigures,"RipleyH_PPSpecific.pdf"), 'ContentType','vector', 'BackgroundColor','none', ...
        'Width',width, 'Height',height, 'Units','centimeters', 'Padding','figure')
end

%% Visualise average scoring patient amongst homogeneous patients (#15)
[sorted_csv_cellTypes,~] = sort(csv_cellTypes);

idxOptimised = minimisationResults.Function2.Both.idxOptimised;
% Keep only homogeneous samples, but stay in sorted order
sortedHomogeneousPatients = sortIdx(ismember(sortIdx,HomogeneousPatients));
np = sortedHomogeneousPatients(round(numel(sortedHomogeneousPatients)./2));
px_name = string(patientInfo.csv_ID(np));

orginalcsv = readmatrix(fullfile(folders.originalCSV,px_name+".csv"));
reps = [3,9,2,4,8,1]; % Randomly generated integers between 1 and 10
options = struct();
options.R = patientInfo.Radius(np);
options.shape="square";
height = 1.88;
width  = 1.88;
options.height = height;
options.width = width;
options.showTicks = false;

%%% Visualise data
figData = figure;
plotConfigurationFct(orginalcsv,sorted_csv_cellTypes,radii,options)

%%% Visualise synthetic data
for ii=1:6
    figData = figure(ii+1);
    synthetic_data = readmatrix(fullfile(folders.SHAPPEOptimisedConfig,px_name,"rep"+num2str(reps(ii))+".csv"));
    plotConfigurationFct(synthetic_data,sorted_csv_cellTypes,radii,options)
end

% Invalid celltypes
valid_cellTypes = false(1,numel(csv_cellTypes));
for nc=1:numel(csv_cellTypes)
    nt = sum(orginalcsv(:,end)==csv_cellTypes(nc));
    if(nt>=minCellCount)
        valid_cellTypes(nc) = true;
    end
end

% Visualise exposure matrix
exportFigureExposure = false;
exportFigureExposureLevels = false;
showColorbar = false;
if(showColorbar)
    figureSuffix = "_colorbar";
else
    figureSuffix = "";
end
statsIMC = neighborhoodExposureStats(orginalcsv(:,[1 2 end]), csv_cellTypes, 50);
Sigma_IMC = statsIMC.e;
Sigma_IMC(~valid_cellTypes,:) = nan;
Sigma_IMC(:,~valid_cellTypes) = nan;
Pi_IMC = statsIMC.mu;
Pi_IMC(~valid_cellTypes,:) = nan;
Pi_IMC(:,~valid_cellTypes) = nan;

Sigma_Synthetic = zeros(size(Sigma_IMC));
Pi_Synthetic = zeros(size(Pi_IMC));

Sigma_Random = zeros(size(Sigma_IMC));
Pi_Random = zeros(size(Pi_IMC));

for rr=1:10
    synthetic_data = readmatrix(fullfile(folders.SHAPPEOptimisedConfig,px_name,"rep"+num2str(rr)+".csv"));
    statsSynthetic = neighborhoodExposureStats(synthetic_data(:,[1 2 end]), csv_cellTypes, 50);
    Sigma_Synthetic = Sigma_Synthetic+statsSynthetic.e;
    Pi_Synthetic = Pi_Synthetic+statsSynthetic.mu;

    random_data = readmatrix(fullfile(folders.ECSRConfigurations,px_name,"rep"+num2str(rr)+".csv"));
    statsRandom =neighborhoodExposureStats(random_data(:,[1 2 end]), csv_cellTypes, 50);
    Sigma_Random = Sigma_Random+statsRandom.e;
    Pi_Random = Pi_Random+statsRandom.mu;
end
Sigma_Synthetic = Sigma_Synthetic./10;
Sigma_Synthetic(~valid_cellTypes,:) = nan;
Sigma_Synthetic(:,~valid_cellTypes) = nan;
Pi_Synthetic = Pi_Synthetic./10;
Pi_Synthetic(~valid_cellTypes,:) = nan;
Pi_Synthetic(:,~valid_cellTypes) = nan;

Sigma_Random = Sigma_Random./10;
Sigma_Random(~valid_cellTypes,:) = nan;
Sigma_Random(:,~valid_cellTypes) = nan;
Pi_Random = Pi_Random./10;
Pi_Random(~valid_cellTypes,:) = nan;
Pi_Random(:,~valid_cellTypes) = nan;

clim_E = [min([min(Sigma_IMC) min(Sigma_Synthetic) min(Sigma_Random)]) ...
    max([min(Sigma_IMC) max(Sigma_Synthetic) max(Sigma_Random)])];
%%% Green shades
c1 = hex2rgb('#F7F5EA');
c2 = hex2rgb('#D6E2C2');
c3 = hex2rgb('#A5C494');
c4 = hex2rgb('#6E9C7B');
c5 = hex2rgb('#41705E');
c6 = hex2rgb('#233F38');
colours = [c1; c2; c3; c4; c5; c6];
breakpoints = clim_E(1)+[0:5].*(clim_E(end)-clim_E(1))./5;
coloursSigma = interp1( breakpoints, colours, linspace(clim_E(1),clim_E(end),256) );

figE = figure;
tiledlayout(1,3,'Padding','Compact')
nexttile
imagesc(Sigma_IMC)
axis square
colormap(coloursSigma)
clim(clim_E)
xticklabels([])
yticklabels([])
set(gca,'FontSize',axis_fontsize,'FontName',axis_fontname,...
    'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
nexttile
imagesc(Sigma_Synthetic)
axis square
colormap(coloursSigma)
clim(clim_E)
xticklabels([])
yticklabels([])
set(gca,'FontSize',axis_fontsize,'FontName',axis_fontname,...
    'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
nexttile
imagesc(Sigma_Random)
axis square
colormap(coloursSigma)
clim(clim_E)
xticklabels([])
yticklabels([])
set(gca,'FontSize',axis_fontsize,'FontName',axis_fontname,...
    'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
if(showColorbar)
    c = colorbar;
    c.Layout.Tile = 'south';
end
height = 2.26;
width  = 7.32;
set(figE, 'Units',  'Centimeters')
set(figE, 'Position',  [0, 0,width,height])
set(figE,'Units','Inches');
pos = get(figE,'Position');
set(figE,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
printFig(figE)
if(exportFigureExposure)
    exportgraphics(figE, fullfile(folders.ClusteringFigures,"exposureHeatMaps"+figureSuffix+".pdf"), 'ContentType','vector', 'BackgroundColor','none', ...
        'Width',width, 'Height',height, 'Units','centimeters', 'Padding','figure')
end

clim_Pi = [min([min(Pi_IMC) min(Pi_Synthetic) min(Pi_Random)]) ...
    max([min(Pi_IMC) max(Pi_Synthetic) max(Pi_Random)])];
breakpoints = clim_Pi(1)+[0:5].*(clim_Pi(end)-clim_Pi(1))./5;
coloursPi = interp1( breakpoints, colours, linspace(clim_Pi(1),clim_Pi(end),256) );
figMu = figure;
tiledlayout(1,3,'Padding','Compact')
nexttile
imagesc(Pi_IMC)
axis square
colormap(coloursPi)
clim(clim_Pi)
xticklabels([])
yticklabels([])
set(gca,'FontSize',axis_fontsize,'FontName',axis_fontname,...
    'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
nexttile
imagesc(Pi_Synthetic)
axis square
colormap(coloursPi)
clim(clim_Pi)
xticklabels([])
yticklabels([])
set(gca,'FontSize',axis_fontsize,'FontName',axis_fontname,...
    'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
nexttile
imagesc(Pi_Random)
axis square
colormap(coloursPi)
clim(clim_Pi)
xticklabels([])
yticklabels([])
set(gca,'FontSize',axis_fontsize,'FontName',axis_fontname,...
    'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
if(showColorbar)
    c = colorbar;
    c.Layout.Tile = 'south';
end
height = 2.26;
width  = 7.32;
set(figMu, 'Units',  'Centimeters')
set(figMu, 'Position',  [0, 0,width,height])
set(figMu,'Units','Inches');
pos = get(figMu,'Position');
set(figMu,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
printFig(figMu)
if(exportFigureExposureLevels)
    exportgraphics(figMu, fullfile(folders.ClusteringFigures,"MuHeatMaps"+figureSuffix+".pdf"), 'ContentType','vector', 'BackgroundColor','none', ...
        'Width',width, 'Height',height, 'Units','centimeters', 'Padding','figure')
end

%% Visualise exposure results
exportFigureBarGraph = true;
%%% Figure 1: Raw H(20), H(40) and H(60)
scoreIndex = [idxSigma,idxPi];
scoreNames = ["\Delta e_{\Sigma}/e_{\Sigma}^{(ECSR)}","\Delta e_{\Pi}/e_{\Pi}^{(ECSR)}"];
figNames = ["Sigma","Pi"];
for ii=1:length(scoreIndex)
    figExposureBarGraph = figure;
    ax = gca;
    allaxes = [ax];
    hold on
    yR = resultsECSR(scoreIndex(ii),sortIdx);
    yH = resultsOptimised(scoreIndex(ii),sortIdx);
    den = yR;
    den(den==0) = 1;
    y = (yR-yH)./den;
    plotIQR(y,figExposureBarGraph,gray);
    b = bar(y,'FaceColor','flat');
    for np=1:Npatients
        if(ismember(sortIdx(np),HomogeneousPatients))
            col = col_homogeneous;
        else
            col = col_inhomogeneous;
        end
        b.CData(np,:) = col;
    end
    yline(0)
    ylabel(scoreNames(ii))
    xticks([])
    yl = ylim;
    ylims(ii,:) = yl;
    ylim(yl)
    xlim([0 Npatients+1])
    set(allaxes,'FontSize',axis_fontsize,'FontName',axis_fontname,...
        'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
    height = 3.23;
    width  = 6.81;
    set(figExposureBarGraph, 'Units',  'Centimeters')
    set(figExposureBarGraph, 'Position',  [0, 0,width,height])
    set(figExposureBarGraph,'Units','Inches');
    pos = get(figExposureBarGraph,'Position');
    set(figExposureBarGraph,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    printFig(figExposureBarGraph)
    if(exportFigureBarGraph)
        exportgraphics(figExposureBarGraph, fullfile(folders.ExposureFigures,figNames(ii)+".pdf"),...
            'ContentType','vector', 'BackgroundColor','none', ...
            'Width',width, 'Height',height, 'Units','centimeters', 'Padding','figure')
    end
end

for ii=1:length(scoreIndex)
    figExposureScatter = figure;
    ax = gca;
    allaxes = [ax];
    yR = resultsECSR(scoreIndex(ii),sortIdx);
    yH = resultsOptimised(scoreIndex(ii),sortIdx);
    hold on
    % Vertical line relying all the random and SHAPPE values
    for i = 1:Npatients
        plot([i, i], [yR(i), yH(i)],'Color', "#3B3B3B", 'LineWidth', 1);
    end
    scatter(1:Npatients,yR,7,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.6,'MarkerEdgeColor',purple)
    scatter(1:Npatients,yH,7,'MarkerFaceColor',orange,'MarkerFaceAlpha',0.6,'MarkerEdgeColor',orange)
    % scatter(sortIdxy(1:2),yR(sortIdxy(1:2)),7,'MarkerFaceColor','black','MarkerFaceAlpha',0.6,'MarkerEdgeColor','black')
    % scatter(sortIdxy(1:2),yR(sortIdxy(1:2)),7,'MarkerFaceColor','black','MarkerFaceAlpha',0.6,'MarkerEdgeColor','black')
    xticks([])
    yl = ylim;
    ylims(ii,:) = yl;
    ylim(yl)
    xlim([0 Npatients+1])
    set(allaxes,'FontSize',axis_fontsize,'FontName',axis_fontname,...
        'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
    height = 3.23;
    width  = 6.81;
    set(figExposureScatter, 'Units',  'Centimeters')
    set(figExposureScatter, 'Position',  [0, 0,width,height])
    set(figExposureScatter,'Units','Inches');
    pos = get(figExposureScatter,'Position');
    set(figExposureScatter,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    printFig(figExposureScatter)
    if(exportFigureBarGraph)
        exportgraphics(figExposureScatter, fullfile(folders.ExposureFigures,figNames(ii)+"_scatter.pdf"),...
            'ContentType','vector', 'BackgroundColor','none', ...
            'Width',width, 'Height',height, 'Units','centimeters', 'Padding','figure')
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DENSITY RESULT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Dissimilarities in Density Maps
exportFigureDM = false;
figDM = figure;
hold on
yR = resultsECSR(idxDMDiff,sortIdx);
yH = resultsOptimised(idxDMDiff,sortIdx);
den = yR;
den(den==0) = 1;
y = (yR-yH)./den;
ciStat = plotIQR(y,figDM,gray);
b = bar(y,'FaceColor','flat','BarWidth',1);
for np=1:Npatients
    if(ismember(sortIdx(np),HomogeneousPatients))
        col = col_homogeneous;
    else
        col = col_inhomogeneous;
    end
    b.CData(np,:) = col;
end
yline(0)
xticks([])
xlabel("IMC Point Patterns")
ylabel("\Delta e_D / e_D^{(ECSR)}")
set(gca,'FontSize',axis_fontsize,'FontName',axis_fontname,...
    'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
height = 3.66;
width  = 13.64;
set(figDM, 'Units',  'Centimeters')
set(figDM, 'Position',  [0, 0,width,height])
set(figDM,'Units','Inches');
pos = get(figDM,'Position');
set(figDM,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
printFig(figDM)
if(exportFigureDM)
    exportgraphics(figDM, fullfile(folders.DensityFigures,"change_eD_SHAPPE.pdf"), 'ContentType','vector', 'BackgroundColor','none', ...
        'Width',width, 'Height',height, 'Units','centimeters', 'Padding','figure')
end

%% Inhomogeneity replication
load(fullfile("inhomogeneity_replicates.mat"))
inhStatus = inhomogeneityIdentification(1,:)>0 & inhomogeneityIdentification(2,:)>0;
inhomogeneity_replicates(:,1) = inhStatus';
scoreInhomogeneityIentification = sum(abs(inhomogeneity_replicates(:,2:end)-inhomogeneity_replicates(:,1))==0,"all")./numel(inhomogeneity_replicates(:,2:end));
exportInhReplicate = false;
%%% Bar graph version
percentInhomogeneityIdentification = 100.*sum(inhomogeneity_replicates(:,2:end),2)./(size(inhomogeneity_replicates,2)-1);
percentInhomogeneityIdentification_inhomogeneous = percentInhomogeneityIdentification(InhomogeneousPatients);
percentInhomogeneityIdentification_homogeneous = percentInhomogeneityIdentification(HomogeneousPatients);

figInhomogeneityReproduction_bargraph = figure;
b = bar([sort(percentInhomogeneityIdentification_inhomogeneous);sort(100-percentInhomogeneityIdentification_homogeneous,'ascend')],...
    'FaceColor','flat','BarWidth',1);
for np=1:Npatients
    if(np<=numel(percentInhomogeneityIdentification_inhomogeneous))
        b.CData(np,:) = col_inhomogeneous;
    else
        b.CData(np,:) = col_homogeneous;
    end
end
yline(50,'LineWidth',plot_linewidth)
xticks(1:Npatients)
xticklabels([])
ylabel("Homogeneity"+newline+"Status Accuracy"+newline+"(%)")
set(gca,'FontSize',axis_fontsize,'FontName',axis_fontname,...
    'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
height = 2.54;
width  = 14.55;
set(figInhomogeneityReproduction_bargraph, 'Units',  'Centimeters')
set(figInhomogeneityReproduction_bargraph, 'Position',  [0, 0,width,height])
set(figInhomogeneityReproduction_bargraph,'Units','Inches');
pos = get(figInhomogeneityReproduction_bargraph,'Position');
set(figInhomogeneityReproduction_bargraph,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
printFig(figInhomogeneityReproduction_bargraph)
if(exportInhReplicate)
    exportgraphics(figInhomogeneityReproduction_bargraph, fullfile(folders.DensityFigures,"inhomogeneityReproduction_barGraph.pdf"), 'ContentType','vector', 'BackgroundColor','none', ...
            'Width',width, 'Height',height, 'Units','centimeters', 'Padding','figure')
end
%% Visualise inhomogeneous patienT
[sorted_csv_cellTypes,sortIdx_celltypes] = sort(csv_cellTypes);
np = 1;
px_name = string(patientInfo.csv_ID(np));
folder_currentFigures = fullfile(folders.DensityFigures,"patient"+string(np));
if(~isfolder(folder_currentFigures))
    mkdir(folder_currentFigures)
end
orginalcsv = readmatrix(fullfile(folders.originalCSV,px_name+".csv"));
inh = ismember(np,InhomogeneousPatients);

reps = randperm(10,6);
options = struct();
options.R = patientInfo.Radius(np);
options.shape="square";
height = 1.88;
width  = 1.88;
options.height = height;
options.width = width;
options.showTicks = false;

%%% Visualise data
figData = figure(1);
plotConfigurationFct(orginalcsv,sorted_csv_cellTypes,radii,options)
printFig(figData)

%%% Visualise synthetic data - Optimised
for ii=1:6
    figData = figure(ii+1);
    synthetic_data = readmatrix(fullfile(folders.SHAPPEOptimisedConfig,px_name,"rep"+num2str(reps(ii))+".csv"));
    plotConfigurationFct(synthetic_data,sorted_csv_cellTypes,radii,options)
end

%%% Visualise synthetic data - Simplified
for ii=1:6
    figData = figure(ii+1+6);
    synthetic_data = readmatrix(fullfile(folders.SHAPPESimplifiedConfig,px_name,"rep"+num2str(reps(ii))+".csv"));
    plotConfigurationFct(synthetic_data,sorted_csv_cellTypes,radii,options)
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALL METRICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

heightHeatmaps = 2.92;
widthHeatmaps = 11.8;

heightCI_IQR = 4.39;
widthCI_IQR = 7.04;

pos_random = 1;
pos_optimised = 3;
pos_add = 5;
pos_optimised_PPSpecific = 7;
pos_larger = 9;
positions = [pos_random pos_optimised pos_add pos_optimised_PPSpecific pos_larger];
a = -0.25;
b = 0.25;

exportMetrics = true;

eNN = [resultsECSR(idxNNDiff,sortIdx);...
    resultsOptimised(idxNNDiff,sortIdx);...
    resultsSimplified(idxNNDiff,sortIdx);...
    resultsPPSpecific(idxNNDiff,sortIdx);...
    resultsDomainIncrease(idxNNDiff,sortIdx)];
eH20 = [resultsECSR(idxH20_rawDiff,sortIdx);...
    resultsOptimised(idxH20_rawDiff,sortIdx);...
    resultsSimplified(idxH20_rawDiff,sortIdx);...
    resultsPPSpecific(idxH20_rawDiff,sortIdx);...
    resultsDomainIncrease(idxH20_rawDiff,sortIdx)];
eH40 = [resultsECSR(idxH40_rawDiff,sortIdx);...
    resultsOptimised(idxH40_rawDiff,sortIdx);...
    resultsSimplified(idxH40_rawDiff,sortIdx);...
    resultsPPSpecific(idxH40_rawDiff,sortIdx);...
    resultsDomainIncrease(idxH40_rawDiff,sortIdx)];
eH60 = [resultsECSR(idxH60_rawDiff,sortIdx);...
    resultsOptimised(idxH60_rawDiff,sortIdx);...
    resultsSimplified(idxH60_rawDiff,sortIdx);...
    resultsPPSpecific(idxH60_rawDiff,sortIdx);...
    resultsDomainIncrease(idxH60_rawDiff,sortIdx)];
eD = [resultsECSR(idxDMDiff,sortIdx);...
    resultsOptimised(idxDMDiff,sortIdx);...
    resultsSimplified(idxDMDiff,sortIdx);...
    resultsPPSpecific(idxDMDiff,sortIdx);...
    resultsDomainIncrease(idxDMDiff,sortIdx)];
eSigma = [resultsECSR(idxSigma,sortIdx);...
    resultsOptimised(idxSigma,sortIdx);...
    resultsSimplified(idxSigma,sortIdx);...
    resultsPPSpecific(idxSigma,sortIdx);...
    resultsDomainIncrease(idxSigma,sortIdx)];
ePi = [resultsECSR(idxPi,sortIdx);...
    resultsOptimised(idxPi,sortIdx);...
    resultsSimplified(idxPi,sortIdx);...
    resultsPPSpecific(idxPi,sortIdx);...
    resultsDomainIncrease(idxPi,sortIdx)];

allMetrics = cat(3,eD,eNN,eH20,eH40,eH60,eSigma,ePi);
figureNames = ["e_D","e_NN","e_H(20)","e_H(40)","e_H(60)","e_Sigma","e_Pi"];

for m=1:size(allMetrics,3)
    metricValues = reshape(allMetrics(:,:,m),size(allMetrics,1:2));
    %%%%% A. Heatmap
    fig = figure;
    ax = axes('Parent',fig,'Position',[0.02 0.28 0.96 0.68]);
    imagesc(ax,metricValues)
    axis(ax,'tight')
    set(ax,'DataAspectRatioMode','auto','PlotBoxAspectRatioMode','auto')
    colormap(ax,flip(seq_colorful_continuous(30:end,:)))
    set(ax,'XTick',[],'YTick',1:5,'YTickLabel',"",'Box','off','Fontsize',7,'FontName',axis_fontname)
    cb = colorbar(ax,'southoutside');
    cb.Position = [0.1 0.13 0.8 0.05];
    set(fig, 'Units',  'Centimeters')
    set(fig, 'Position',  [0, 0,widthHeatmaps,heightHeatmaps])
    set(fig,'Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    printFig(fig)
    if(exportMetrics)
        exportgraphics(fig,fullfile(folders.MetricsFigures,figureNames(m)+".pdf"),'Resolution',600,...
            'BackgroundColor', 'none', 'ContentType', 'vector','Width',widthHeatmaps, 'Height',heightHeatmaps,...
            'Units','centimeters', 'Padding','figure')
    end

    %%%%% B. 95% CI of median
    posVector = zeros([numel(positions) size(metricValues,2)]);
    fig_95CI = figure;
    hold on
    for p=1:numel(positions)
        data = metricValues(p,:);
        ci = bootci(1000, {@median, data});
        % Scatter of metric
        tempPos = rand(1, size(eNN,2));
        posVector(p,:) = tempPos;
        scatter(positions(p).*ones(size(metricValues(1,:)))+a+(b-a).*posVector(p,:),...
            metricValues(p,:),8,"MarkerFaceColor",gray,"MarkerFaceAlpha",0.6,"MarkerEdgeColor",gray);
        % Median
        plot([positions(p)+a positions(p)+b],[median(data) median(data)],'Color',dark_red,'LineWidth',plot_linewidth)
        % CI
        h=fill([positions(p)+0.5.*a,positions(p)+0.5.*b,...
            positions(p)+0.5.*b,positions(p)+0.5.*a],[ci(1),ci(1),ci(2),ci(2)],dark_red,'EdgeColor',dark_red);
        h.FaceAlpha=0.6;
    end
    xticks(positions)
    xticklabels([])
    set(fig_95CI, 'Units',  'Centimeters')
    set(fig_95CI, 'Position',  [0, 0,widthCI_IQR,heightCI_IQR])
    set(fig_95CI,'Units','Inches');
    pos = get(fig_95CI,'Position');
    set(fig_95CI,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    printFig(fig_95CI)
    if(exportMetrics)
        exportgraphics(fig_95CI,fullfile(folders.MetricsFigures,figureNames(m)+"_95CI.pdf"),'Resolution',600,...
            'BackgroundColor', 'none', 'ContentType', 'vector','Width',widthCI_IQR, 'Height',heightCI_IQR,...
            'Units','centimeters', 'Padding','figure')
    end
    %%%%% C. IQR
    height = 4.39;
    width = 7.04;
    fig_IQR = figure;
    hold on
    for p=1:numel(positions)
        data = metricValues(p,:);
        IQR = quantile(data,[0.25 0.75]);
        % Scatter of metric
        scatter(positions(p).*ones(size(metricValues(1,:)))+a+(b-a).*posVector(p,:),...
            metricValues(p,:),8,"MarkerFaceColor",gray,"MarkerFaceAlpha",0.6,"MarkerEdgeColor",gray);
        % Median
        plot([positions(p)+a positions(p)+b],[median(data) median(data)],'Color',dark_blue,'LineWidth',plot_linewidth)
        % CI
        h=fill([positions(p)+0.5.*a,positions(p)+0.5.*b,...
            positions(p)+0.5.*b,positions(p)+0.5.*a],[IQR(1),IQR(1),IQR(2),IQR(2)],dark_blue,'EdgeColor',dark_blue);
        h.FaceAlpha=0.6;
    end
    xticks(positions)
    xticklabels([])
    set(fig_IQR, 'Units',  'Centimeters')
    set(fig_IQR, 'Position',  [0, 0,widthCI_IQR,heightCI_IQR])
    set(fig_IQR,'Units','Inches');
    pos = get(fig_IQR,'Position');
    set(fig_IQR,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    printFig(fig_IQR)
    if(exportMetrics)
        exportgraphics(fig_IQR,fullfile(folders.MetricsFigures,figureNames(m)+"_IQR.pdf"),'Resolution',600,...
            'BackgroundColor', 'none', 'ContentType', 'vector','Width',widthCI_IQR, 'Height',heightCI_IQR,...
            'Units','centimeters', 'Padding','figure')
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% HYPERPARAMETER IMPORTANCE %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pareto Front
% paramComb = [K(:), S(:), RW(:), N(:), WNoCorr(:), NAlpha(:)];   % columns: k, s, pw, N

scoresdata = [F2_avg,F3_avg,F4_avg];
idxPareto = minimisationResults.Function2.Both.idxParetoMembership;
idxNonPareto = setdiff(1:nComb,idxPareto);
colPareto = repmat(purples(3,:),numel(idxPareto),1);
colNonPareto = repmat(light_blues(4,:),numel(idxNonPareto),1);

exportParetoFront = false;
figPareto = figure;
s = scatter3(F2_avg(idxNonPareto),F3_avg(idxNonPareto),F4_avg(idxNonPareto),0.6*marker_size,colNonPareto,'MarkerFaceColor','flat',...
    'MarkerFaceAlpha',0.6);
hold on
s2 = scatter3(F2_avg(idxPareto),F3_avg(idxPareto),F4_avg(idxPareto),marker_size,colPareto,'MarkerFaceColor','flat',...
    'MarkerFaceAlpha',0.6);
xlabel("e_{H}(20)")
ylabel("e_{H}(40)")
zlabel("e_{H}(60)")
width = 5.18;
height = 4.06;
set(gca,'FontSize',axis_fontsize,'FontName',axis_fontname,...
    'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
set(figPareto, 'Units',  'Centimeters')
set(figPareto, 'Position',  [0, 0,width,height])
set(figPareto,'Units','Inches');
posFig = get(figPareto,'Position');
set(figPareto,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[posFig(3), posFig(4)])
printFig(figPareto)
if(exportParetoFront)
    exportgraphics(figPareto,fullfile(folders.HyperparameterOptimisationFigures,"ParetoFrontier.pdf"),'Resolution',600,...
        'BackgroundColor', 'none', 'ContentType', 'vector','Width',width, 'Height',height,...
        'Units','centimeters', 'Padding','figure')
end

[R,pval] = corrcoef(scoresdata);

%% Importance

hyperparametersImportance_filename = "hyperparameterImportance.csv";
hyperparametersImportance = readmatrix(hyperparametersImportance_filename);

exportFigureHyperparameterImportance = false;
figImportance = figure;
b = bar(hyperparametersImportance,'BarWidth',1);
b(1).FaceColor = light_blues(4,:);
b(2).FaceColor = purples(1,:);
b(3).FaceColor = purples(3,:);
b(4).FaceColor = purples(4,:);
width = 9.14;
height = 4.06;
xticklabels([])
set(gca,'FontSize',axis_fontsize,'FontName',axis_fontname,...
    'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
set(figImportance, 'Units',  'Centimeters')
set(figImportance, 'Position',  [0, 0,width,height])
set(figImportance,'Units','Inches');
posFig = get(figImportance,'Position');
set(figImportance,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[posFig(3), posFig(4)])
printFig(figImportance)
if(exportFigureHyperparameterImportance)
    exportgraphics(figImportance,fullfile(folders.HyperparameterOptimisationFigures,"importance.pdf"),'Resolution',600,...
        'BackgroundColor', 'none', 'ContentType', 'vector','Width',width, 'Height',height,...
        'Units','centimeters', 'Padding','figure')
end

%%% For the close up
% xl=xlim;
% xlim([xl(1) 3.5])
% ylim([0 0.003])
% yticks(ylim)
% newWidth = 3.51;
% newHeight = 1.12;
% set(figImportance, 'Units',  'Centimeters')
% set(figImportance, 'Position',  [0, 0,newWidth,newHeight])
% set(figImportance,'Units','Inches');
% posFig = get(figImportance,'Position');
% set(figImportance,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[posFig(3), posFig(4)])
% printFig(figImportance)
% if(exportFigureHyperparameterImportance)
%     exportgraphics(figImportance,fullfile(folders.HyperparameterOptimisationFigures,"importance_closeup.pdf"),'Resolution',600,...
%         'BackgroundColor', 'none', 'ContentType', 'vector','Width',width, 'Height',height,...
%         'Units','centimeters', 'Padding','figure')
% end

%% Perform a chi square test for parameters
% Copula - phi - rw - L - wNoCorr - alpha
pvalsChiGOF = ones(1,size(paramComb,2));
paramsNames = ["Copula","phi","rw","L","wNoCorr","alpha"];
for pp=1:size(paramComb,2)
    disp(paramsNames(pp))
    % Get unique values
    vals = unique(paramComb(:,pp));

    % We need to treat alpha differently because alpha is undefined for L=1
    if pp==6
        optimisedValues = optimisedSets_PPSpecific(pp,optimisedSets_PPSpecific(4,:)~=1);
    else
        optimisedValues = optimisedSets_PPSpecific(pp,:);
    end

    % Calculate expected counts: How many parameter combinations have this
    % value times the number of samples
    expCounts = zeros(1,numel(vals));
    realCounts = zeros(1,numel(vals));
    for v=1:numel(vals)
        if(pp==6)
            expCounts(1,v) = numel(optimisedValues).*sum(paramComb(paramComb(:,4)~=1,pp)==vals(v))./size(paramComb(paramComb(:,4)~=1,pp),1);
        else
            expCounts(1,v) = numel(optimisedValues).*sum(paramComb(:,pp)==vals(v))./nComb;
        end
        realCounts(1,v) = sum(optimisedValues==vals(v));
    end
    disp("Expected counts:")
    disp(expCounts)
    disp("Real counts:")
    disp(realCounts)
    chi2stat = sum((realCounts - expCounts).^2 ./ expCounts);
    p = 1 - chi2cdf(chi2stat, numel(vals) - 1);
    resid = (realCounts - expCounts) ./ sqrt(expCounts);
    disp("Standardised residuals:"); disp(resid)
    pvalsChiGOF(1,pp) = p;
end

% Correction for multiple testing
[adjusted_p_values] = mafdr(pvalsChiGOF,'BHFDR',true);

%% Visualise frequency of hyperparameter values
exportHyperparametersDistributions = false;
colors = [light_blues(5,:);light_blues(2,:);light_blues(4,:)];
kVals  = [1 2]; % kernel
sVals  = [2 5 10]; % step limit
rwVals = [1 5 10]; % position weight
NVals  = [1 5 10]; % Number of neighbors
wNoCorrVals = [0, 5, 10]; % 
neighborAlphaVals = [0, 0.5, 1];
% [K(:), S(:), RW(:), N(:), WNoCorr(:), NAlpha(:)];   % columns: k, s, pw, N

% Kernel - step limit - repulsion weight - L - wNoCorr - alpha
% Significant before correction : Kernel - step limit - repulsion weight -
% wNoCorr
% Significant after correction : repulsion weight - wNoCorr
chiSquareResults = zeros(1,6);

space_inbetween = 2;
pos = space_inbetween;
ticks = [];
labels_figure = [];
labels = ["All","Hom.","Inh."];

figHyperparameters = figure;
hold on
%%% Kernel
vals = zeros(3,3);
optisedKernels = optimisedSets_PPSpecific(1,:);
countsHom = zeros(numel(kVals),1);
countsInh = zeros(numel(kVals),1);
for i=1:numel(kVals)
    vals(1,i) = sum(optisedKernels==kVals(i))./Npatients;
    countsHom(i,1) = sum(optisedKernels(HomogeneousPatients)==kVals(i));
    countsInh(i,1) = sum(optisedKernels(InhomogeneousPatients)==kVals(i));
    vals(2,i) = countsHom(i,1)./numel(HomogeneousPatients);
    vals(3,i) = countsInh(i,1)./numel(InhomogeneousPatients);
end
[~,pval,~] = FisherExactTest([countsHom,countsInh]);
chiSquareResults(1,1) = pval;
for i=1:3
    b = bar(pos,100.*vals(i,:),'stacked','BarWidth',0.9);
    ticks = [ticks,pos];
    labels_figure = [labels_figure,labels(i)]; 
    for c=1:length(b)
        b(c).FaceColor = colors(c,:);
    end
    pos = pos+1;
end
pos = pos+space_inbetween;

%%% Step limit
vals = zeros(3,3);
optisedValues = optimisedSets_PPSpecific(2,:);
countsHom = zeros(numel(sVals),1);
countsInh = zeros(numel(sVals),1);
for i=1:numel(sVals)
    vals(1,i) = sum(optisedValues==sVals(i))./Npatients;
    countsHom(i,1) = sum(optisedValues(HomogeneousPatients)==sVals(i));
    countsInh(i,1) = sum(optisedValues(InhomogeneousPatients)==sVals(i));
    vals(2,i) = countsHom(i,1)./numel(HomogeneousPatients);
    vals(3,i) = countsInh(i,1)./numel(InhomogeneousPatients);
end
[~,pval,~] = FisherExactTest([countsHom,countsInh]);
chiSquareResults(1,2) = pval;
for i=1:3
    b = bar(pos,100.*vals(i,:),'stacked','BarWidth',0.9);
    ticks = [ticks,pos];
    labels_figure = [labels_figure,labels(i)]; 
    for c=1:length(b)
        b(c).FaceColor = colors(c,:);
    end
    pos = pos+1;
end
pos = pos+space_inbetween;

%%% Repulsion weight
vals = zeros(3,3);
optisedValues = optimisedSets_PPSpecific(3,:);
countsHom = zeros(numel(rwVals),1);
countsInh = zeros(numel(rwVals),1);
for i=1:numel(rwVals)
    vals(1,i) = sum(optisedValues==rwVals(i))./Npatients;
    countsHom(i,1) = sum(optisedValues(HomogeneousPatients)==rwVals(i));
    countsInh(i,1) = sum(optisedValues(InhomogeneousPatients)==rwVals(i));
    vals(2,i) = countsHom(i,1)./numel(HomogeneousPatients);
    vals(3,i) = countsInh(i,1)./numel(InhomogeneousPatients);
end
[~,pval,~] = FisherExactTest([countsHom,countsInh]);
chiSquareResults(1,3) = pval;
for i=1:3
    b = bar(pos,100.*vals(i,:),'stacked','BarWidth',0.9);
    ticks = [ticks,pos];
    labels_figure = [labels_figure,labels(i)]; 
    for c=1:length(b)
        b(c).FaceColor = colors(c,:);
    end
    pos = pos+1;
end
pos = pos+space_inbetween;

%%% Number of neighbors
vals = zeros(3,3);
optisedValues = optimisedSets_PPSpecific(4,:);
countsHom = zeros(numel(NVals),1);
countsInh = zeros(numel(NVals),1);
for i=1:numel(NVals)
    vals(1,i) = sum(optisedValues==NVals(i))./Npatients;
    countsHom(i,1) = sum(optisedValues(HomogeneousPatients)==NVals(i));
    countsInh(i,1) = sum(optisedValues(InhomogeneousPatients)==NVals(i));
    vals(2,i) = countsHom(i,1)./numel(HomogeneousPatients);
    vals(3,i) = countsInh(i,1)./numel(InhomogeneousPatients);
end
[~,pval,~] = FisherExactTest([countsHom,countsInh]);
chiSquareResults(1,4) = pval;
for i=1:3
    b = bar(pos,100.*vals(i,:),'stacked','BarWidth',0.9);
    ticks = [ticks,pos];
    labels_figure = [labels_figure,labels(i)]; 
    for c=1:length(b)
        b(c).FaceColor = colors(c,:);
    end
    pos = pos+1;
end
pos = pos+space_inbetween;

%%% Weight no correlation
vals = zeros(3,3);
optisedValues = optimisedSets_PPSpecific(5,:);
countsHom = zeros(numel(wNoCorrVals),1);
countsInh = zeros(numel(wNoCorrVals),1);
for i=1:numel(wNoCorrVals)
    vals(1,i) = sum(optisedValues==wNoCorrVals(i))./Npatients;
    countsHom(i,1) = sum(optisedValues(HomogeneousPatients)==wNoCorrVals(i));
    countsInh(i,1) = sum(optisedValues(InhomogeneousPatients)==wNoCorrVals(i));
    vals(2,i) = countsHom(i,1)./numel(HomogeneousPatients);
    vals(3,i) = countsInh(i,1)./numel(InhomogeneousPatients);
end
[~,pval,~] = FisherExactTest([countsHom,countsInh]);
chiSquareResults(1,5) = pval;
for i=1:3
    b = bar(pos,100.*vals(i,:),'stacked','BarWidth',0.9);
    ticks = [ticks,pos];
    labels_figure = [labels_figure,labels(i)]; 
    for c=1:length(b)
        b(c).FaceColor = colors(c,:);
    end
    pos = pos+1;
end
pos = pos+space_inbetween;

%%% Weight of neighbors (alpha)
vals = zeros(3,3);
optisedValues = optimisedSets_PPSpecific(6,:);
countsHom = zeros(numel(neighborAlphaVals),1);
countsInh = zeros(numel(neighborAlphaVals),1);
for i=1:numel(neighborAlphaVals)
    vals(1,i) = sum(optisedValues==neighborAlphaVals(i))./Npatients;
    countsHom(i,1) = sum(optisedValues(HomogeneousPatients)==neighborAlphaVals(i));
    countsInh(i,1) = sum(optisedValues(InhomogeneousPatients)==neighborAlphaVals(i));
    vals(2,i) = countsHom(i,1)./numel(HomogeneousPatients);
    vals(3,i) = countsInh(i,1)./numel(InhomogeneousPatients);
end
[~,pval,~] = FisherExactTest([countsHom,countsInh]);
chiSquareResults(1,6) = pval;
for i=1:3
    b = bar(pos,100.*vals(i,:),'stacked','BarWidth',0.9);
    ticks = [ticks,pos];
    labels_figure = [labels_figure,labels(i)]; 
    for c=1:length(b)
        b(c).FaceColor = colors(c,:);
    end
    pos = pos+1;
end
xlim([0 pos+space_inbetween-1])
width = 14.66;
height = 3.66;
xticks(ticks)
xticklabels(labels_figure)
ylabel("Representation in"+newline+"Optimised Sets (%)")
set(gca,'FontSize',axis_fontsize,'FontName',axis_fontname,...
    'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
set(figHyperparameters, 'Units',  'Centimeters')
set(figHyperparameters, 'Position',  [0, 0,width,height])
set(figHyperparameters,'Units','Inches');
posFig = get(figHyperparameters,'Position');
set(figHyperparameters,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[posFig(3), posFig(4)])
printFig(figHyperparameters)
if(exportHyperparametersDistributions)
    exportgraphics(figHyperparameters,fullfile(folders.HyperparameterOptimisationFigures,"Distributions.pdf"),'Resolution',600,...
        'BackgroundColor', 'none', 'ContentType', 'vector','Width',width, 'Height',height,...
        'Units','centimeters', 'Padding','figure')
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Execution times %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

executionTimes = readmatrix(fullfile(folders.SHAPPEOptimisedConfig,"executionTimes.csv"));
% Read convergence information
convergence1 = readmatrix(fullfile(folders.SHAPPEOptimisedConfig,"convergence1.csv"));
convergence2 = readmatrix(fullfile(folders.SHAPPEOptimisedConfig,"convergence2.csv"));
convergence3 = readmatrix(fullfile(folders.SHAPPEOptimisedConfig,"convergence3.csv"));
convergence4 = readmatrix(fullfile(folders.SHAPPEOptimisedConfig,"convergence4.csv"));
%%executionTimes = zeros(nComb,4); % min - max - mean - med
convergences = zeros(Ncelltypes,4);
for nc=1:Ncelltypes
    convergences(nc,1) = sum(convergence1(~isnan(convergence1(:,nc)),nc));
    convergences(nc,2) = sum(convergence2(~isnan(convergence2(:,nc)),nc));
    convergences(nc,3) = sum(convergence3(~isnan(convergence3(:,nc)),nc));
    convergences(nc,4) = sum(convergence4(~isnan(convergence4(:,nc)),nc));
end % End of for nc=1:Ncelltypes
convergences = convergences(:,1:3);
convergences_all = sum(convergences,1)./sum(convergences(:));
convergences_celltypes = convergences./sum(convergences,2);

CriteriaCols = [orange;dark_blue;dark_red];

exportConvergenceFigure = true;
width = 5.52;
height = 3.25;
figConvergence = figure;
b = bar(convergences_all,'FaceColor','flat','BarWidth',1,'FaceAlpha',0.75);
for bb=1:size(convergences_all,2)
    b.CData(bb,:) = CriteriaCols(bb,:);
end
xticklabels([]);
set(gca,'FontSize',axis_fontsize,'FontName',axis_fontname,...
    'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
set(figConvergence, 'Units',  'Centimeters')
set(figConvergence, 'Position',  [0, 0,width,height])
set(figConvergence,'Units','Inches');
posFig = get(figConvergence,'Position');
set(figConvergence,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[posFig(3), posFig(4)])
printFig(figConvergence)
if(exportConvergenceFigure)
    exportgraphics(figConvergence,fullfile(folders.ExecutionConvergenceFigures,"Convergence.pdf"),'Resolution',600,...
        'BackgroundColor', 'none', 'ContentType', 'vector','Width',width, 'Height',height,...
        'Units','centimeters', 'Padding','figure')
end

exportConvergenceCellTypesFigure = true;
width = 7.41;
height = 3.25;
figConvergenceCellTypes = figure;
b = bar(convergences_celltypes,'stacked','FaceColor','flat','BarWidth',1);
for bb=1:numel(b)
    b(bb).FaceColor = CriteriaCols(bb,:);
    b(bb).FaceAlpha = 0.75;
end
xticklabels([]);
set(gca,'FontSize',axis_fontsize,'FontName',axis_fontname,...
    'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
set(figConvergenceCellTypes, 'Units',  'Centimeters')
set(figConvergenceCellTypes, 'Position',  [0, 0,width,height])
set(figConvergenceCellTypes,'Units','Inches');
posFig = get(figConvergenceCellTypes,'Position');
set(figConvergenceCellTypes,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[posFig(3), posFig(4)])
printFig(figConvergenceCellTypes)
if(exportConvergenceCellTypesFigure)
    exportgraphics(figConvergenceCellTypes,fullfile(folders.ExecutionConvergenceFigures,"ConvergenceCellTypes.pdf"),'Resolution',600,...
        'BackgroundColor', 'none', 'ContentType', 'vector','Width',width, 'Height',height,...
        'Units','centimeters', 'Padding','figure')
end

%%% Close up
figConvergenceCellTypesCloseup = figure;
b = bar(convergences_celltypes,'stacked','FaceColor','flat','BarWidth',1);
for bb=1:numel(b)
    b(bb).FaceColor = CriteriaCols(bb,:);
    b(bb).FaceAlpha = 0.75;
end
xticklabels([]);
ylim([0.75 1])
yticks([0.9 1])
set(gca,'FontSize',axis_fontsize,'FontName',axis_fontname,...
    'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
set(figConvergenceCellTypesCloseup, 'Units',  'Centimeters')
set(figConvergenceCellTypesCloseup, 'Position',  [0, 0,width,height])
set(figConvergenceCellTypesCloseup,'Units','Inches');
posFig = get(figConvergenceCellTypesCloseup,'Position');
set(figConvergenceCellTypesCloseup,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[posFig(3), posFig(4)])
printFig(figConvergenceCellTypesCloseup)
if(exportConvergenceCellTypesFigure)
    exportgraphics(figConvergenceCellTypesCloseup,fullfile(folders.ExecutionConvergenceFigures,"ConvergenceCellTypes_closeup.pdf"),'Resolution',600,...
        'BackgroundColor', 'none', 'ContentType', 'vector','Width',width, 'Height',height,...
        'Units','centimeters', 'Padding','figure')
end

%% % Plot execution times
exportExecutionFigure = true;

% Rebuild the order since it was shuffled
allJobs = [repelem((1:Npatients)',Nreplicates), repmat((1:Nreplicates)',Npatients,1)];
rng(123456,'twister');
allJobs = allJobs(randperm(size(allJobs,1)),:);
[~,order] = sortrows(allJobs);
executionTimes_sorted = executionTimes(order);

% One row per patient, one column per replicate
timesByPatient = reshape(executionTimes_sorted, Nreplicates, Npatients)';
a = -0.2;
b = 0.2;

figExecutionTimes = figure;
width  = 13.67;
height = 3.25;
hold on
% Plot a line for the median
plot([1 Npatients],[median(executionTimes) median(executionTimes)]./60,'LineWidth',plot_linewidth,...
    'Color',pink)
% Plot IQR
iqr = quantile(executionTimes./60,[0.25 0.75]);
h=fill([0,Npatients+1,Npatients+1,0],[iqr(1),iqr(1),iqr(2),iqr(2)],pink,'EdgeColor','none');
h.FaceAlpha=0.6;
for p=1:Npatients
    data = timesByPatient(p,:)./60;
    tempPos = rand(1, Nreplicates);
    scatter(p.*ones(size(data))+a+(b-a).*tempPos,...
        data,8,"MarkerFaceColor",gray,"MarkerFaceAlpha",0.6,"MarkerEdgeColor",gray);
end
xlim([0 Npatients+1])
xticks(1:Npatients)
xticklabels([])
xlabel("")
ylabel("")
set(figExecutionTimes, 'Units',  'Centimeters')
set(figExecutionTimes, 'Position',  [0, 0,width,height])
set(figExecutionTimes,'Units','Inches');
pos = get(figExecutionTimes,'Position');
set(figExecutionTimes,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
ax = gca;
ax.Units = 'normalized';
ti = ax.TightInset + 0.01;
ax.Position = [ti(1), ti(2), 1-ti(1)-ti(3), 1-ti(2)-ti(4)];
printFig(figExecutionTimes)
if(exportExecutionFigure)
    exportgraphics(figExecutionTimes,fullfile(folders.ExecutionConvergenceFigures,"executionTimes.pdf"),'Resolution',600,...
        'BackgroundColor', 'none', 'ContentType', 'vector','Width',width, 'Height',height,...
        'Units','centimeters', 'Padding','tight')
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Calculate the proportion of inside the CI bounds %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load confidence intervals
LBs = loadVariable("RipleyLBs_markBootstrap.mat");
UBs = loadVariable("RipleyUBs_markBootstrap.mat");

insideBoundsSHAPPE = zeros(Npatients,numel(rL));
insideBoundsimplifiedSHAPPE = zeros(Npatients,numel(rL));
insideBoundsRandom = zeros(Npatients,numel(rL));

warning('off', 'all')

for np=1:Npatients

    disp(np)
    patientname = string(patientInfo.csv_ID(np));
    R = patientInfo.Radius(np);

    % Load original data
    IMC = readmatrix(fullfile(folders.originalCSV,patientname+".csv"));

    %%% Visualise Ripley's H function
    OriginalData.folder = folders.originalCSV;
    OriginalData.celltypes = csv_cellTypes;
    OriginalData.color = dark_blue;
    OriginalData.name = "IMC Data";
    OriginalData.shape = "square";
    
    SHAPPE.folder = folders.SHAPPEOptimisedConfig;
    SHAPPE.celltypes = csv_cellTypes;
    SHAPPE.color = orange;
    SHAPPE.name = "SHAPPE";
    SHAPPE.shape = "square";

    SHAPPESimplified.folder = folders.SHAPPESimplifiedConfig;
    SHAPPESimplified.celltypes = csv_cellTypes;
    SHAPPESimplified.color = pink;
    SHAPPESimplified.name = "SIMPLIFIED SHAPPE";
    SHAPPESimplified.shape = "square";
    
    ECSR.folder = folders.ECSRConfigurations;
    ECSR.celltypes = csv_cellTypes;
    ECSR.color = purple;
    ECSR.name = "Random";
    ECSR.shape = "square";
       
    [crossFig,rawScores] = crossCorrelationPlots(patientname,...
        'FolderPatientsData',pwd,...
        'OriginalData',OriginalData,...
        'Method1',SHAPPE,...
        'Method2',SHAPPESimplified,...
        'Method3',ECSR,...
        'PatientInfoTable',patientInfo,...
        'CellTypeNames',["GBM","CD8+","Stroma","M0","M1","M2"],...
        'ReturnScores',false, ...
        'HidePlot',false,...
        'RipleyDistances',rL,...
        'compareMethods',false,...
        'homogeneous',ismember(np,HomogeneousPatients),...
        'Nreps',10,...
        'markerSize',7,...
        'returnInsideScores',true);
    close(crossFig)

    num_validComparisons = num_validComparisons+rawScores.numValidComparisons;
    insideBoundsSHAPPE(np,:) = rawScores.M1./rawScores.numValidComparisons;
    insideBoundsimplifiedSHAPPE(np,:) = rawScores.M2./rawScores.numValidComparisons;
    insideBoundsRandom(np,:) = rawScores.M3./rawScores.numValidComparisons;

end

save('insideBoundsSHAPPE','insideBoundsSHAPPE')
save('insideBoundsimplifiedSHAPPE','insideBoundsimplifiedSHAPPE')
save('insideBoundsRandom','insideBoundsRandom')
%%
load("insideBoundsSHAPPE.mat")
load('insideBoundsimplifiedSHAPPE.mat')
load('insideBoundsRandom.mat')

sortScoresBounds = mean(insideBoundsSHAPPE,2);
[~,sortIdxBounds] = sort(sortScoresBounds);

col_random = purple;
col_optimised = orange;
col_simplified = dark_blue;
col_PPSpecific = pink;
col_larger = light_blue;
marker_size_bounds = 1.*marker_size;

exportFigureVAS = true;
figBounds20 = figure;
yline(0.5,'LineWidth',plot_linewidth,'Color',gray)
hold on
scatter(1:Npatients,insideBoundsSHAPPE(sortIdx,1),marker_size_bounds,'Marker','o',...
    'MarkerFaceColor',col_optimised,'MarkerEdgeColor',col_optimised,'MarkerFaceAlpha',0.6)
scatter(1:Npatients,insideBoundsRandom(sortIdx,1),marker_size_bounds,'Marker','o',...
    'MarkerFaceColor',col_random,'MarkerEdgeColor',col_random,'MarkerFaceAlpha',0.6)
width = 4.52;
height = 2.94;
xticks(1:Npatients)
xticklabels([])
xlim([0 Npatients+1])
ylabel("VAS(20\mum)")
xlabel("IMC Point Pattern")
set(gca,'FontSize',axis_fontsize,'FontName',axis_fontname,...
    'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
set(figBounds20, 'Units',  'Centimeters')
set(figBounds20, 'Position',  [0, 0,width,height])
set(figBounds20,'Units','Inches');
posFig = get(figBounds20,'Position');
set(figBounds20,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[posFig(3), posFig(4)])
printFig(figBounds20)
if(exportFigureVAS)
    exportgraphics(figBounds20,fullfile(folders.InsideBoundsFigures,"20.pdf"),'Resolution',600,...
        'BackgroundColor', 'none', 'ContentType', 'vector','Width',width, 'Height',height,...
        'Units','centimeters', 'Padding','figure')
end

figBounds40 = figure;
yline(0.5,'LineWidth',plot_linewidth,'Color',gray)
hold on
scatter(1:Npatients,insideBoundsSHAPPE(sortIdx,2),marker_size_bounds,'Marker','o',...
    'MarkerFaceColor',col_optimised,'MarkerEdgeColor',col_optimised,'MarkerFaceAlpha',0.6)
scatter(1:Npatients,insideBoundsRandom(sortIdx,2),marker_size_bounds,'Marker','o',...
    'MarkerFaceColor',col_random,'MarkerEdgeColor',col_random,'MarkerFaceAlpha',0.6)
width = 4.52;
height = 2.94;
xticks(1:Npatients)
xticklabels([])
xlim([0 Npatients+1])
ylabel("VAS(40\mum)")
xlabel("IMC Point Pattern")
set(gca,'FontSize',axis_fontsize,'FontName',axis_fontname,...
    'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
set(figBounds40, 'Units',  'Centimeters')
set(figBounds40, 'Position',  [0, 0,width,height])
set(figBounds40,'Units','Inches');
posFig = get(figBounds40,'Position');
set(figBounds40,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[posFig(3), posFig(4)])
printFig(figBounds40)
if(exportFigureVAS)
    exportgraphics(figBounds40,fullfile(folders.InsideBoundsFigures,"40.pdf"),'Resolution',600,...
        'BackgroundColor', 'none', 'ContentType', 'vector','Width',width, 'Height',height,...
        'Units','centimeters', 'Padding','figure')
end

figBounds60 = figure;
yline(0.5,'LineWidth',plot_linewidth,'Color',gray)
hold on
scatter(1:Npatients,insideBoundsSHAPPE(sortIdx,3),marker_size_bounds,'Marker','o',...
    'MarkerFaceColor',col_optimised,'MarkerEdgeColor',col_optimised,'MarkerFaceAlpha',0.6)
scatter(1:Npatients,insideBoundsRandom(sortIdx,3),marker_size_bounds,'Marker','o',...
    'MarkerFaceColor',col_random,'MarkerEdgeColor',col_random,'MarkerFaceAlpha',0.6)
width = 4.52;
height = 2.94;
xticks(1:Npatients)
xticklabels([])
xlim([0 Npatients+1])
ylabel("VAS(60\mum)")
xlabel("IMC Point Pattern")
set(gca,'FontSize',axis_fontsize,'FontName',axis_fontname,...
    'XColor',box_color,'YColor',box_color,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1)
set(figBounds60, 'Units',  'Centimeters')
set(figBounds60, 'Position',  [0, 0,width,height])
set(figBounds60,'Units','Inches');
posFig = get(figBounds60,'Position');
set(figBounds60,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[posFig(3), posFig(4)])
printFig(figBounds60)
if(exportFigureVAS)
    exportgraphics(figBounds60,fullfile(folders.InsideBoundsFigures,"60.pdf"),'Resolution',600,...
        'BackgroundColor', 'none', 'ContentType', 'vector','Width',width, 'Height',height,...
        'Units','centimeters', 'Padding','figure')
end




