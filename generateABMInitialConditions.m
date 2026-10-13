close all
clearvars

%% Retrieve paths to folder

folders.code = pwd;
folders.originalIMC = "/Volumes/SeagateRecherche/ProjectMechanisticReproductionPointPatterns/patientsData/OriginalCSV";
folders.ECSRReconstructions_IMC = "/Volumes/SeagateRecherche/ProjectMechanisticReproductionPointPatterns/"+...
    "Reconstructions_EnlargedDomains/ECSR_Reconstructions_IMCtypes";
folders.SHAPPEReconstructions_IMC = "/Volumes/SeagateRecherche/ProjectMechanisticReproductionPointPatterns/"+...
    "Reconstructions_EnlargedDomains/Reconstructions_IMCtypes";
folders.RepulsiveRegionsSHAPPE = "/Volumes/SeagateRecherche/ProjectMechanisticReproductionPointPatterns/"+...
    "Reconstructions_EnlargedDomains/AugmentedVoidInfoCircles";

folders.ECSRReconstructions_ABM = "/Volumes/SeagateRecherche/ProjectMechanisticReproductionPointPatterns/"+...
    "Reconstructions_EnlargedDomains/ECSR_Reconstructions_ABMtypes";
folders.SHAPPEReconstructions_ABM = "/Volumes/SeagateRecherche/ProjectMechanisticReproductionPointPatterns/"+...
    "Reconstructions_EnlargedDomains/Reconstructions_ABMtypes";
folders.CSRReconstructions_IMC = "/Volumes/SeagateRecherche/ProjectMechanisticReproductionPointPatterns/"+...
    "Reconstructions_EnlargedDomains/CSR_Reconstructions_IMCtypes";
folders.CSRReconstructions_ABM = "/Volumes/SeagateRecherche/ProjectMechanisticReproductionPointPatterns/"+...
    "Reconstructions_EnlargedDomains/CSR_Reconstructions_ABMtypes";
folders.RepulsiveRegionsABM = "/Volumes/SeagateRecherche/ProjectMechanisticReproductionPointPatterns/"+...
    "Reconstructions_EnlargedDomains/RepulsiveRegions";

% Make sure all folders are created
fields = string(fieldnames(folders));
for ff=1:length(fieldnames(folders))
    f = folders.(fields(ff));
    if(~isfolder(f))
        mkdir(f)
    end
end

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

Nreps = 10;
Ncells = 10000;

%% Load file with patients' information
patientInfo_together = readtable(fullfile("IMC_TMEReconstructions",'SOCPatientsInfo.csv'));
Npatients = size(patientInfo_together,1);

%% Load file with inhomogeneity information
inhomogeneityIdentification = loadVariable(fullfile("IMC_TMEReconstructions","inhomogeneityFinalResults_50.mat"));
inhomogeneouSamples = find((inhomogeneityIdentification(1,:)>0) & (inhomogeneityIdentification(2,:)>0));

%% Generate Augmented Configurations

for np=1:Npatients
    patientName = string(patientInfo_together.csv_ID(np));
    inhomogeneityStatus = ismember(np,inhomogeneouSamples);
    NcellsOriginal = patientInfo_together.number_of_cells(np);
    RObs = patientInfo_together.Radius(np);

    % Make sure all folders are created
    if(~isfolder(fullfile(folders.SHAPPEReconstructions_ABM,patientName)))
        mkdir(fullfile(folders.SHAPPEReconstructions_ABM,patientName))
    end
    if(~isfolder(fullfile(folders.ECSRReconstructions_ABM,patientName)))
        mkdir(fullfile(folders.ECSRReconstructions_ABM,patientName))
    end
    if(~isfolder(fullfile(folders.CSRReconstructions_IMC,patientName)))
        mkdir(fullfile(folders.CSRReconstructions_IMC,patientName))
    end
    if(~isfolder(fullfile(folders.CSRReconstructions_ABM,patientName)))
        mkdir(fullfile(folders.CSRReconstructions_ABM,patientName))
    end

    % Map SHAPPE Reconstructions to ABM types
    for nr=1:Nreps
        config = readmatrix(fullfile(folders.SHAPPEReconstructions_IMC,patientName,"rep"+string(nr)+".csv"));
        new_config = mapIMC_to_ABM(config);
        writematrix(new_config,fullfile(folders.SHAPPEReconstructions_ABM,patientName,"rep"+string(nr)+".csv"))
    end
    % Map ECSR Reconstructions to ABM types
    for nr=1:Nreps
        config = readmatrix(fullfile(folders.ECSRReconstructions_IMC,patientName,"rep"+string(nr)+".csv"));
        new_config = mapIMC_to_ABM(config);
        writematrix(new_config,fullfile(folders.ECSRReconstructions_ABM,patientName,"rep"+string(nr)+".csv"))
        % If sample was homogeneous, ECSR and CSR reconstructions are the same
        if(~inhomogeneityStatus)
            writematrix(config,fullfile(folders.CSRReconstructions_IMC,patientName,"rep"+string(nr)+".csv"))
            writematrix(new_config,fullfile(folders.CSRReconstructions_ABM,patientName,"rep"+string(nr)+".csv"))
        end
    end
    % If sample was inhomogeneous, we need to generate configurations for
    % CSR
    if(inhomogeneityStatus)
        initialProportions = zeros(1,Ncelltypes);
        for t=1:Ncelltypes
            initialProportions(t) = patientInfo_together.("proportion_"+celltypes_name(t))(np);
        end
        for nr=1:Nreps
            [CSR_IMC,~,~,~,~] = ECSRGenerator(Ncells,initialProportions,Ncells./NcellsOriginal,RObs,"circle","square",...
                0.8202,csv_cellTypes,true);
            CSR_ABM = mapIMC_to_ABM(CSR_IMC);
            writematrix(CSR_IMC,fullfile(folders.CSRReconstructions_IMC,patientName,"rep"+string(nr)+".csv"))
            writematrix(CSR_ABM,fullfile(folders.CSRReconstructions_ABM,patientName,"rep"+string(nr)+".csv"))
        end
    end
end


%% AUGMENTED VOIDS: Generating voids for very large random configurations (Attracting and Repulsing)
thresholdInformation = struct();
thresholdInformation.DensityThreshold = 0;
bL = 40;
max_tries = 1000;
for np=1:Npatients
    patientName =string(patientInfo_together.csv_ID(np));
    inhomogeneityStatus = ismember(np,inhomogeneouSamples);
    NcellsOriginal = patientInfo_together.number_of_cells(np);
    RObs = patientInfo_together.Radius(np);

    % Make sure folder is created
    if(~isfolder(fullfile(folders.RepulsiveRegionsABM,patientName)))
        mkdir(fullfile(folders.RepulsiveRegionsABM,patientName))
    end

    initialProportions = zeros(1,Ncelltypes);
    for t=1:Ncelltypes
        initialProportions(csv_cellTypes(t)) = patientInfo_together.(...
            "proportion_"+celltypes_name(t))(np);
    end

    % Load IMC Data
    originalConfig = readmatrix(fullfile(folders.originalIMC,patientName+".csv"));
    % Perform analysis of zero-density regions
    zero_density_analysis = densityGridAnalysis(originalConfig,RObs,bL,"square",thresholdInformation);
    % Increase domain
    R_A_initial = patientInfo_together.synthetic_radius_circle(np);
    R_A = 2.*R_A_initial;
    areaABM = pi.*(R_A.^2);
    areaIMC = (2.*RObs)^2;
    for nr=1:Nreps
        % Load voids of initial configurations
        voids_initial = readmatrix(fullfile(folders.RepulsiveRegionsSHAPPE,patientName,"rep"+string(nr)+".csv"));
        if(size(voids_initial,1)>0)
            %%% How many cold sections?
            NcoldIMC = zero_density_analysis.num_patches;
            coldAreas = zero_density_analysis.AreaDistribution;
            Ncold = round(NcoldIMC.*(areaABM/areaIMC));
            % Generate voids
            voids = zeros(Ncold,3);
            voids(1:size(voids_initial,1),:) = voids_initial;
            placeVoids = true;
            v = size(voids_initial,1)+1;
            current_try = 1;
            while(v<=Ncold)
                % Sample effective radius
                if(v<=NcoldIMC)
                    area = coldAreas(v);
                else
                    area = coldAreas(randi(length(coldAreas)));
                end
                r = sqrt(area./pi);
                maxRadius = R_A-r;
                minRadius = R_A_initial+r;
                while(minRadius>=maxRadius)
                    area = coldAreas(randi(length(coldAreas)));
                    r = sqrt(area./pi);
                    maxRadius = R_A-r;
                    minRadius = R_A_initial+r;
                end
                % Sample center
                thetaPos = 2*pi*rand();
                rPos = sqrt(rand()*(minRadius^2-maxRadius^2)+maxRadius^2);
                x = rPos*cos(thetaPos);
                y = rPos*sin(thetaPos);
                if(v==1)
                    voids(v,:) = [x,y,r];
                    v = v+1;
                else
                    % Check if overlapping with already placed voids
                    ok = true;
                    for vv=1:(v-1)
                        x_other = voids(vv,1);
                        y_other = voids(vv,2);
                        r_other = voids(vv,3);
                        d = norm([x,y]-[x_other,y_other]);
                        if(d<=(r_other+r+20))
                            ok = false;
                            break
                        end % End of if(d<=(r_other+r))
                    end % End of for vv=1:(v-1)
                    if(ok)
                        voids(v,:) = [x,y,r];
                        v = v+1;
                        current_try = 1;
                    else
                        current_try = current_try+1;
                        if(current_try>max_tries)
                            current_try = 1;
                            v = size(voids_initial,1);
                        end % End of if(current_try>max_tries)
                    end % End of if(ok)
                end % End of if(v==1)
            end % End of while(v<Ncold)
        else
            voids = voids_initial;
        end
        writematrix(voids,fullfile(folders.RepulsiveRegionsABM,patientName,"rep"+string(nr)+".csv"))
    end

end


%% Calculate probability of being inside a void
voidProbs = zeros(1,Npatients);
voidAverageDensities = zeros(1,Npatients);
patients = 1:Npatients;
for np=patients
    if(ismember(np,inhomogeneouSamples) && np~=2)
        px_name = string(patientInfo_together.csv_ID(np));
        probs = zeros(1,Nreps);
        avgDensities = zeros(1,Nreps);
        px_configurationFolder = fullfile(folders.SHAPPEReconstructions_IMC,px_name);
        px_voidsFolder = fullfile(folders.RepulsiveRegionsSHAPPE,px_name);
        for nr=1:Nreps
            % Load configuration
            data = readmatrix(fullfile(px_configurationFolder,"rep"+string(nr)+".csv"));
            Ncells = size(data,1);
            % Load voids
            voids = readmatrix(fullfile(px_voidsFolder,"rep"+string(nr)+".csv"));
            % Calculate how many cells are inside a void
            cellcounts = 0;
            for c=1:Ncells
                x = data(c,1);
                y = data(c,2);
                d2 = sqrt((voids(:,1)-x).^2+(voids(:,2)-y).^2)-voids(:,3);
                if(any(d2<0))
                    cellcounts = cellcounts+1;
                end
            end
            probs(1,nr) = cellcounts./Ncells;

            % Calculate density inside voids
            voids_density = zeros(1,size(voids,1));
            for nv=1:size(voids,1)
                % Calculate how many cells are inside the void
                d2 = sqrt((voids(nv,1)-data(:,1)).^2+(voids(nv,2)-data(:,2)).^2)-voids(nv,3);
                cellcounts = sum(d2<0);
                voids_density(1,nv) = cellcounts./(pi.*(voids(nv,3)).^2);
            end
            if(size(voids,1)>0)
                inside_density = mean(voids_density);
            end
            avgDensities(1,nr) = inside_density;
        end
        voidProbs(1,np) = mean(probs);
        voidAverageDensities(1,np) = mean(avgDensities);
    end
end

writematrix(voidProbs,fullfile("/Volumes/SeagateRecherche/ProjectMechanisticReproductionPointPatterns/Reconstructions_EnlargedDomains/voidProbs_New.csv"))
writematrix(voidAverageDensities,fullfile("/Volumes/SeagateRecherche/ProjectMechanisticReproductionPointPatterns/Reconstructions_EnlargedDomains/voidProbsAreaUnit_New.csv"))

%% Calculate cell density inside voids
configuration_folder = fullfile(pwd,"algorithmConfigurations_N10000",...
    "HookeConfigurations_Gaussian_stepLimit_10_repulsionWeight_1_Neigh_10_weightRandom_0_neighAlpha_0/");
voids_folder = fullfile(pwd,"patientsData",...
    "AugmentedVoidInfoCircles/");
%"AugmentedVoidInfoCircles_3foldInitialDomain/");
voidProbs = zeros(1,Npatients);
patients = 1:Npatients;
Nreps = 10;
for np=patients
    if(~ismember(np,homogeneousPatients) && np~=2)
        px_name = string(patientInfo_together.csv_ID(np));
        probs = zeros(1,Nreps);
        px_configurationFolder = fullfile(configuration_folder,px_name);
        px_voidsFolder = fullfile(voids_folder,px_name);
        for nr=1:Nreps
            % Load configuration
            data = readmatrix(fullfile(px_configurationFolder,"rep"+string(nr)+".csv"));
            Ncells = size(data,1);
            % Load voids
            voids = readmatrix(fullfile(px_voidsFolder,"rep"+string(nr)+".csv"));
            voids_density = zeros(1,size(voids,1));
            for nv=1:size(voids,1)
                % Calculate how many cells are inside the void
                d2 = sqrt((voids(nv,1)-data(:,1)).^2+(voids(nv,2)-data(:,2)).^2)-voids(nv,3);
                cellcounts = sum(d2<0);
                voids_density(1,nv) = cellcounts./(pi.*(voids(nv,3)).^2);
            end
            if(size(voids,1)>0)
                inside_density = mean(voids_density);
            end
            probs(1,nr) = inside_density;
        end
        voidProbs(1,np) = mean(probs);
    end
end

writematrix(voidProbs,fullfile(pwd,"patientsData","voidProbsAreaUnit.csv"))