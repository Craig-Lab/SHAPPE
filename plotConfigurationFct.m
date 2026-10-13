function plotConfigurationFct(coord,cellTypes,radius,options)

    if(nargin<3)
        marker_size = 10;
        xL = [min(coord(:,1))-0.1.*abs(min(coord(:,1))) max(coord(:,1))+0.1*abs(max(coord(:,1)))];
        yL = [min(coord(:,2))-0.1.*abs(min(coord(:,2))) max(coord(:,2))+0.1*abs(max(coord(:,2)))];
        R = max(max(abs(xL)),max(abs(yL)));
        set_colors;
        colors = [oranges(3,:);pinks(3,:);purples(3,:);dark_blues(3,:);light_blues(3,:);dark_reds(3,:)];
        plotLgd = false;
        radius = [];
        shape="square";
    elseif(nargin<4)
        marker_size = 10;
        xL = [min(coord(:,1))-0.1.*abs(min(coord(:,1))) max(coord(:,1))+0.1*abs(max(coord(:,1)))];
        yL = [min(coord(:,2))-0.1.*abs(min(coord(:,2))) max(coord(:,2))+0.1*abs(max(coord(:,2)))];
        set_colors;
        colors = [oranges(3,:);pinks(3,:);purples(3,:);dark_blues(3,:);light_blues(3,:);dark_reds(3,:)];
        plotLgd = false;
        R = max(max(abs(xL)),max(abs(yL)));
        shape="square";
    else
        if(isfield(options,'markersize'))
            marker_size = options.markersize;
        else
            marker_size = 36;
        end
        
        if(isfield(options,'R'))
            R = options.R;
        else
            R = max(max(abs(xL)),max(abs(yL)));
        end

        if(isfield(options,'xlim'))
            xL = options.xlim;
        else
            xL = [-ceil(R) ceil(R)];
        end
    
        if(isfield(options,'ylim'))
            yL = options.ylim;
        else
            yL = [-ceil(R) ceil(R)];
        end
        
        if(isfield(options,'colors'))
            colors = options.colors;
        else
            set_colors;
            colors = [oranges(3,:);pinks(3,:);purples(3,:);dark_blues(3,:);light_blues(3,:);dark_reds(3,:)];
        end
        
        if(isfield(options,'plotLegend'))
            plotLgd = options.plotLegend;
            cellTypesLegend = options.cellTypesLegend;
        else
            plotLgd = false;
        end

        if(isfield(options,'shape'))
            shape = options.shape;
        else
            shape = "square";
        end
        if(isfield(options,"width"))
            width = options.width;
        else
            width = [];
        end
        if(isfield(options,"height"))
            height = options.height;
        else
            height = [];
        end
        if(isfield(options,"showTicks"))
            showTicks = options.showTicks;
        else
            showTicks = true;
        end
    end

    hold on
    % Plot a bounding box
    patch([xL(1) xL(end) xL(end) xL(1)], ...
      [yL(1) yL(1) yL(end) yL(end)], ...
      'w', ...
      'FaceAlpha', 0, ...
      'EdgeColor', 'none');
    if(strcmp(shape,"square"))
        patch([-R R R -R -R], ...
          [-R -R R R -R],'w','EdgeColor', 'none');
        plot([-R R R -R -R],[-R -R R R -R],'Color','black','LineWidth',0.75)
    elseif(strcmp(shape,"circle"))
        theta = linspace(0,2.*pi.*300);
        patch(R.*cos(theta), ...
          R.*sin(theta),'w','EdgeColor', 'none');
        plot(R.*cos(theta),R.*sin(theta),'Color','black','LineWidth',0.75)
    else
        errMsg = "Only shape=square or shape=circle are supported for the domain boundaries.";
        error(errMsg)
    end
    for t=1:length(cellTypes)
        h{t} = scatter(coord(coord(:,end)==cellTypes(t),1),coord(coord(:,end)==cellTypes(t),2),marker_size,...
            'MarkerEdgeColor',colors(t,:),'MarkerFaceColor',colors(t,:),'MarkerFaceAlpha',0.6,'LineWidth',0.2);
    end
    if(plotLgd)
        legend(cellTypesLegend)
    end
    
    axis square
    if(isempty(xL))
        xlim([-ceil(R) ceil(R)])
        ylim([-ceil(R) ceil(R)])
        xticks([-ceil(R) ceil(R)])
        yticks([-ceil(R) ceil(R)])
    else
        xlim(xL)
        ylim(yL)
    end
    if(~showTicks)
        xticks([])
        yticks([])
    end

    if(~isempty(height))
        set(gcf, 'Units',  'Centimeters')
        set(gcf, 'Position',  [0, 0,width,height])
        set(gcf,'Units','Inches');
        pos = get(gcf,'Position');
        set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    end

    if(~isempty(radius))
        % Obtain the axes size (in axpos) in Points
        currentunits = get(gca,'Units');
        set(gca, 'Units', 'Points');
        axpos = get(gca,'Position');
        set(gca, 'Units', currentunits);
        for t=1:length(cellTypes)
            desiredDiameter = radius(t).*2;
            markerWidth = desiredDiameter/diff(xlim)*axpos(3); % Calculate Marker width in points
            set(h{t}, 'SizeData', markerWidth^2)
        end
    end

    printFig(gcf)
    
end