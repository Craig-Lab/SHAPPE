function varargout = crossCorrelationPlots(csv_name,varargin)

    p = inputParser;
    
    validFolder = @(x) (ischar(x) || isstring(x)) && isfolder(x);
    validFolderOrEmpty = @(x) isempty(x) || ((ischar(x) || isstring(x)) && isfolder(x));
    
    addRequired(p,'csv_name',@(x) ischar(x) || isstring(x));
    
    addParameter(p,'FolderPatientsData',"",validFolder);
    addParameter(p,'OriginalData',"");
    addParameter(p,'Method1',struct([]));
    addParameter(p,'Method2',struct([]));
    addParameter(p,'Method3',struct([]));
    addParameter(p,'PatientInfoTable',[]);
    addParameter(p,'CellTypeNames',[]);
    addParameter(p,'ReturnScores',false);
    addParameter(p,'HidePlot',false);
    addParameter(p,'RipleyDistances',[20]);
    addParameter(p,'Nreps',5);
    addParameter(p,'minCellCount',15);
    addParameter(p,'compareMethods',false);
    addParameter(p,'homogeneous',true);
    addParameter(p,'figSize',[]);
    addParameter(p,'markerSize',20);
    addParameter(p,'returnInsideScores',false);
    
    parse(p, csv_name, varargin{:});

    FolderPatientsData = p.Results.FolderPatientsData;
    OriginalData = p.Results.OriginalData;
    returnInsideScores = p.Results.returnInsideScores;
    Method1 = p.Results.Method1;
    CellTypeNames = p.Results.CellTypeNames;
    ReturnScores = p.Results.ReturnScores;
    HidePlot = p.Results.HidePlot;
    compareMethods = p.Results.compareMethods;
    homogeneous = p.Results.homogeneous;
    figSize = p.Results.figSize;
    markerSize = p.Results.markerSize;
    plotMethod2 = false;
    if(~isempty(p.Results.Method2))
        plotMethod2 = true;
        Method2 = p.Results.Method2;
    end
    
    plotMethod3 = false;
    if(~isempty(p.Results.Method3))
        plotMethod3 = true;
        Method3 = p.Results.Method3;
    end
    PatientInfoTable = p.Results.PatientInfoTable;
    % Retrieve the number associated to the patient
    np = find(strcmp(PatientInfoTable.csv_ID(:),csv_name));
    RadiusData = PatientInfoTable.Radius(np);
    if(~isfield(Method1,"Radius"))
        RadiusNewConfiguration = PatientInfoTable.Radius(np);
    else
        RadiusNewConfiguration = Method1.Radius;
    end

    % Retrieve bounds
    rL = p.Results.RipleyDistances;
    if ReturnScores || compareMethods
        nMethods = 1 + plotMethod2 + plotMethod3;
        scores = nan(nMethods, length(rL));
    end

    Ncelltypes = length(CellTypeNames);
    Nreps = p.Results.Nreps;
    LBs = load(fullfile(FolderPatientsData,"RipleyLBs_markBootstrap"));
    LBs = LBs.LBRipley_allPatients;
    UBs = load(fullfile(FolderPatientsData,"RipleyUBs_markBootstrap"));
    UBs = UBs.UBRipley_allPatients;
    minCellCount = p.Results.minCellCount;
    LB_CurrentPatient= reshape(LBs(np,:,:,:),[Ncelltypes,Ncelltypes,length(rL)]);
    UB_CurrentPatient = reshape(UBs(np,:,:,:),[Ncelltypes,Ncelltypes,length(rL)]);
    

    RipleyData = nan(Ncelltypes,Ncelltypes,length(rL));
    RipleyMethod1 = nan(Ncelltypes,Ncelltypes,length(rL),Nreps);
    RipleyMethod2 = nan(Ncelltypes,Ncelltypes,length(rL),Nreps);
    RipleyMethod3 = nan(Ncelltypes,Ncelltypes,length(rL),Nreps);

    dataOriginalConfig = readmatrix(fullfile(OriginalData.folder,csv_name+".csv"));
    LdiffscoreMethod1 = 0;
    LdiffscoreMethod2 = 0;
    LdiffscoreMethod3 = 0;
    insideScoreM1 = zeros(1,length(rL));
    numValidComparisons = 0;

    CIwidth = UB_CurrentPatient-LB_CurrentPatient;
    allScoresMethod1 = nan(length(rL),Ncelltypes,Ncelltypes,Nreps);
    if(plotMethod2)
        allScoresMethod2 = nan(length(rL),Ncelltypes,Ncelltypes,Nreps);
        insideScoreM2 = zeros(1,length(rL));;
    end
    if(plotMethod3)
        allScoresMethod3 = nan(length(rL),Ncelltypes,Ncelltypes,Nreps);
        insideScoreM3 = zeros(1,length(rL));;
    end

    for r=1:Nreps
        LdiffscoreMethod1_rep = 0;
        LdiffscoreMethod2_rep = 0;
        LdiffscoreMethod3_rep = 0;
        sumweights = 0;
        dataMethod1 = readmatrix(fullfile(Method1.folder,csv_name,"rep"+string(r)+".csv"));
        if(plotMethod2)
            dataMethod2 = readmatrix(fullfile(Method2.folder,csv_name,"rep"+string(r)+".csv"));
        end % End of if(plotMethod2)
        if(plotMethod3)
            dataMethod3 = readmatrix(fullfile(Method3.folder,csv_name,"rep"+string(r)+".csv"));
        end % End of if(plotMethod3)

        % Loop through cell types
        for nc=1:Ncelltypes
            coordOriginal = dataOriginalConfig(dataOriginalConfig(:,end)==OriginalData.celltypes(nc),1:2);
            coordMethod1 = dataMethod1(dataMethod1(:,end)==Method1.celltypes(nc),1:2);
            if(plotMethod2)
                coordMethod2 = dataMethod2(dataMethod2(:,end)==Method2.celltypes(nc),1:2);
            end % End of if(plotMethod2)
            if(plotMethod3)
                coordMethod3 = dataMethod3(dataMethod3(:,end)==Method3.celltypes(nc),1:2);
            end % End of if(plotMethod2)
            if(size(coordOriginal,1)>=minCellCount)
                for nc2=1:Ncelltypes
                    Nj = length(find(dataOriginalConfig(:,end)==OriginalData.celltypes(nc2)));
                    if(Nj>=minCellCount)
                        sumweights = sumweights+1;
                        numValidComparisons = numValidComparisons+1;
                        if(r==1)
                            if(homogeneous)
                                [L,~] = ripleyK(dataOriginalConfig,OriginalData.celltypes(nc),OriginalData.celltypes(nc2),RadiusData,rL,"H" ,OriginalData.shape);
                            else
                                [L,~] = ripleyKInhom(dataOriginalConfig,OriginalData.celltypes(nc),OriginalData.celltypes(nc2),RadiusData,rL,100,"H" ,OriginalData.shape);
                            end
                            RipleyData(nc,nc2,:) = L;
                        end
                        LB = sqrt(max(0,reshape(LB_CurrentPatient(nc,nc2,:),[1,length(rL)]))./pi)-rL;
                        UB = sqrt(reshape(UB_CurrentPatient(nc,nc2,:),[1,length(rL)])./pi)-rL;
                        if(homogeneous)
                            [LMethod1,~] = ripleyK(dataMethod1,Method1.celltypes(nc),Method1.celltypes(nc2),RadiusNewConfiguration,rL,"H" ,Method1.shape);
                            if(plotMethod2)
                                [LMethod2,~] = ripleyK(dataMethod2,Method2.celltypes(nc),Method2.celltypes(nc2),RadiusNewConfiguration,rL,"H" ,Method2.shape);
                                RipleyMethod2(nc,nc2,:,r) = LMethod2;
                            end
                            if(plotMethod3)
                                [LMethod3,~] = ripleyK(dataMethod3,Method3.celltypes(nc),Method3.celltypes(nc2),RadiusNewConfiguration,rL,"H" ,Method3.shape);
                                RipleyMethod3(nc,nc2,:,r) = LMethod3;
                            end
                        else
                            [LMethod1,~] = ripleyKInhom(dataMethod1,Method1.celltypes(nc),Method1.celltypes(nc2),RadiusNewConfiguration,rL,100,"H" ,Method1.shape);
                            if(plotMethod2)
                                [LMethod2,~] = ripleyKInhom(dataMethod2,Method2.celltypes(nc),Method2.celltypes(nc2),RadiusNewConfiguration,rL,100,"H" ,Method2.shape);
                                RipleyMethod2(nc,nc2,:,r) = LMethod2;
                            end
                            if(plotMethod3)
                                [LMethod3,~] = ripleyKInhom(dataMethod3,Method3.celltypes(nc),Method3.celltypes(nc2),RadiusNewConfiguration,rL,100,"H" ,Method3.shape);
                                RipleyMethod3(nc,nc2,:,r) = LMethod3;
                            end
                        end
                        RipleyMethod1(nc,nc2,:,r) = LMethod1;
                        % scoresMethod1 = zeros(1,length(rL));
                        % scoresMethod1(LMethod1>=LB & LMethod1<=UB) = 1;
                        scoresMethod1 = abs(L-LMethod1);
                        LdiffscoreMethod1_rep = LdiffscoreMethod1_rep+scoresMethod1;

                        inside = zeros(1,length(rL));
                        inside(LMethod1>=LB & LMethod1<=UB) = 1;
                        insideScoreM1 = insideScoreM1+inside;

                        den = UB-LB;
                        den(den==0) = 1;
                        exceedanceMethod1 = (max(0, LB - LMethod1) + max(0, LMethod1 - UB))./den;
                        allScoresMethod1(:,nc,nc2,r) = exceedanceMethod1(:);

                        if(plotMethod2)
                            % scoresMethod2 = zeros(1,length(rL));
                            % scoresMethod2(LMethod2>=LB & LMethod2<=UB) = 1;
                            % insideScoreM2 = insideScoreM2+scoresMethod2;

                            scoresMethod2 = abs(L-LMethod2);
                            LdiffscoreMethod2_rep = LdiffscoreMethod2_rep+scoresMethod2;
    
                            inside = zeros(1,length(rL));
                            inside(LMethod2>=LB & LMethod2<=UB) = 1;
                            insideScoreM2 = insideScoreM2+inside;

                            exceedanceMethod2 = (max(0, LB - LMethod2) + max(0, LMethod2 - UB))./den;
                            allScoresMethod2(:,nc,nc2,r) = exceedanceMethod2(:);
                        end
                        if(plotMethod3)
                            % scoresMethod3 = zeros(1,length(rL));
                            % scoresMethod3(LMethod3>=LB & LMethod3<=UB) = 1;
                            % insideScoreM3 = insideScoreM3+scoresMethod3;

                            scoresMethod3 = abs(L-LMethod3);
                            LdiffscoreMethod3_rep = LdiffscoreMethod3_rep+scoresMethod3;
    
                            inside = zeros(1,length(rL));
                            inside(LMethod3>=LB & LMethod3<=UB) = 1;
                            insideScoreM3 = insideScoreM3+inside;

                            exceedanceMethod3 = (max(0, LB - LMethod3) + max(0, LMethod3 - UB))./den;
                            allScoresMethod3(:,nc,nc2,r) = exceedanceMethod3(:);
                        end
                    end % end of if(Nj>=minCellCount)
                end % End of for nc2=1:Ncelltypes
            end
        end % for nc=1:Ncelltypes
        LdiffscoreMethod1 = LdiffscoreMethod1+LdiffscoreMethod1_rep./sumweights;
        if(plotMethod2)
            LdiffscoreMethod2 = LdiffscoreMethod2+LdiffscoreMethod2_rep./sumweights;
        end
        if(plotMethod3)
            LdiffscoreMethod3 = LdiffscoreMethod3+LdiffscoreMethod3_rep./sumweights;
        end
    end % End of for r=1:Nreps

    LdiffscoreMethod1 = LdiffscoreMethod1./Nreps;
    allScoresMethod1 = mean(allScoresMethod1,4);
    if(ReturnScores)
        scores(1,:) = LdiffscoreMethod1;
    end
    if(plotMethod2)
        LdiffscoreMethod2 = LdiffscoreMethod2./Nreps;
        if(ReturnScores)
            scores(2,:) = LdiffscoreMethod2;
        end
        allScoresMethod2 = mean(allScoresMethod2,4);
    end
    if(plotMethod3)
        LdiffscoreMethod3 = LdiffscoreMethod3./Nreps;
        if(ReturnScores)
            scores(3,:) = LdiffscoreMethod3;
        end
        allScoresMethod3 = mean(allScoresMethod3,4);
    end
    if(compareMethods)
        comparisons = nan(size(scores,1)-1,length(rL),Ncelltypes,Ncelltypes);
        for d=1:length(rL)
            for nc1=1:Ncelltypes
                for nc2=1:Ncelltypes
                    if(~isnan(allScoresMethod1(d,nc1,nc2)))
                        idxComparison = 1;
                        tau = 0.1;
                        e1 = allScoresMethod1(d,nc1,nc2);
                        if(plotMethod2)
                            e2 = allScoresMethod2(d,nc1,nc2);
                            if(e1<(e2-tau))
                                % Win
                                comparisons(idxComparison,d,nc1,nc2) = 1;
                            elseif(abs(e1-e2)<=tau)
                                % Tie
                                comparisons(idxComparison,d,nc1,nc2) = 2;
                            else
                                % Loss
                                comparisons(idxComparison,d,nc1,nc2) = 3;
                            end
                            idxComparison = idxComparison+1;
                        end
                        if(plotMethod3)
                            e3 = allScoresMethod3(d,nc1,nc2);
                            if(e1<(e3-tau))
                                % Win
                                comparisons(idxComparison,d,nc1,nc2) = 1;
                            elseif(abs(e1-e3)<=tau)
                                % Tie
                                comparisons(idxComparison,d,nc1,nc2) = 2;
                            else
                                % Loss
                                comparisons(idxComparison,d,nc1,nc2) = 3;
                            end
                        end
                    end
                end
            end
        end
    end

    if(~HidePlot)
        colBounds = 169.*ones(1,3)./255;
        % Generate plot
        fig = figure;
        tl = tiledlayout(Ncelltypes,Ncelltypes);
        nexttile
        hBounds   = plot(nan,nan,'-','LineWidth',1.5,'Color',colBounds); hold on
        hOriginal = scatter(nan,nan,markerSize,'filled','MarkerFaceColor',OriginalData.color,'MarkerEdgeColor',OriginalData.color,...
            'MarkerFaceAlpha',0.6,'Marker','diamond');
        hM1       = scatter(nan,nan,markerSize,'filled','MarkerFaceColor',Method1.color,'MarkerEdgeColor',Method1.color,'MarkerFaceAlpha',0.6);
        if(plotMethod2)
            hM2 = scatter(nan,nan,markerSize,'filled','MarkerFaceColor',Method2.color,'MarkerEdgeColor',Method2.color,'MarkerFaceAlpha',0.6);
        end
        if(plotMethod3)
            hM3 = scatter(nan,nan,markerSize,'filled','MarkerFaceColor',Method3.color,'MarkerEdgeColor',Method3.color,'MarkerFaceAlpha',0.6);
        end
        for nc1=1:Ncelltypes
            for nc2=1:Ncelltypes
                nexttile((nc1-1).*Ncelltypes+nc2)
                if(nc1==1)
                    title(CellTypeNames(nc2))
                end
                if(nc2==1)
                    ylabel(CellTypeNames(nc1))
                    set(get(gca,'YLabel'),'Rotation',0)
                end
                if(~isnan(RipleyData(nc1,nc2,1)))
                    LB = sqrt(max(0,reshape(LB_CurrentPatient(nc1,nc2,:),[1,length(rL)]))./pi)-rL;
                    UB = sqrt(reshape(UB_CurrentPatient(nc1,nc2,:),[1,length(rL)])./pi)-rL;
                    hold on
                    % Plot Lower and Upper bounds
                    errorbar(rL,(LB+UB)./2,(UB-LB)./2,'Color',colBounds,'LineWidth',1)
                   
                    % Plot Method 1
                    for r=1:Nreps
                        scatter(rL,reshape(RipleyMethod1(nc1,nc2,:,r),[1,length(rL)]),markerSize,...
                            'MarkerFaceColor',Method1.color,'MarkerEdgeColor',Method1.color,'MarkerFaceAlpha',0.6)
                    end
                    % Plot Method 2
                    if(plotMethod2)
                        for r=1:Nreps
                            scatter(rL,reshape(RipleyMethod2(nc1,nc2,:,r),[1,length(rL)]),markerSize,...
                                'MarkerFaceColor',Method2.color,'MarkerEdgeColor',Method2.color,'MarkerFaceAlpha',0.6)
                        end
                    end
                    % PLot Method 3
                    if(plotMethod3)
                        for r=1:Nreps
                            scatter(rL,reshape(RipleyMethod3(nc1,nc2,:,r),[1,length(rL)]),markerSize,...
                                'MarkerFaceColor',Method3.color,'MarkerEdgeColor',Method3.color,'MarkerFaceAlpha',0.6)
                        end
                    end
                    % Plot Data
                    scatter(rL,squeeze(RipleyData(nc1,nc2,:)),2.*markerSize,...
                        'MarkerFaceColor',OriginalData.color,'MarkerEdgeColor',OriginalData.color,'MarkerFaceAlpha',1,'Marker','diamond')
                else
                    text(0.5, 0.5, 'N/A', ...
                        'Units', 'normalized', ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle', ...
                        'FontSize', 14, ...
                        'FontWeight', 'bold');
                end % End of if(~isnan(RipleyData(nc1,nc2,1)))
            end % End of for nc2=1:Ncelltypes
        end % End of for nc1=1:Ncelltypes
        lgdArticles = [hBounds hOriginal hM1];
        lgdNames = ["95% CI" OriginalData.name Method1.name];
        if(plotMethod2)
            lgdArticles = [lgdArticles, hM2];
            lgdNames = [lgdNames, Method2.name];
        end
        if(plotMethod3)
            lgdArticles = [lgdArticles, hM3];
            lgdNames = [lgdNames, Method3.name];
        end
        lgd = legend(lgdArticles,lgdNames,'NumColumns',length(lgdNames));
        lgd.Layout.Tile = 'south';
    end

    varargout{1} = fig;
    k = 2;
    if(ReturnScores)
        varargout{k} = scores;
        k = k+1;
    end
    if(compareMethods)
        varargout{k} = comparisons;
        k = k+1;
    end
    if(returnInsideScores)
        RawScores = struct();
        RawScores.M1 = insideScoreM1;
        if(plotMethod2)
            RawScores.M2 = insideScoreM2;
        end
        if(plotMethod3)
            RawScores.M3 = insideScoreM3;
        end
        RawScores.numValidComparisons = numValidComparisons;
        varargout{k} = RawScores;
    end
end
