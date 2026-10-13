function F = ripleyF(pointPattern,R,shape,distances,varargin)

    p = inputParser;
    addRequired(p,'pointPattern');
    addRequired(p,'R');
    addRequired(p,'shape');
    addRequired(p,'distances');
    addParameter(p,'NTestLocations',10000);

    parse(p, pointPattern, R,shape, distances,varargin{:});
    NTestLocations = p.Results.NTestLocations;

    % Generate random points in the domain
    if(strcmp(shape,"circle"))
        theta = 2*pi*rand(NTestLocations,1);
        r = R*sqrt(rand(NTestLocations,1));
        testPositions = [r.*cos(theta) r.*sin(theta)];
    elseif(strcmp(shape,"square"))
        testPositions = unifrnd(-R,R,NTestLocations,2);
    else
        msg = "Shape must be circle or square.";
        error(msg)
    end
    % Build a kd-tree to search nearest neighbors
    kdtree = KDTreeSearcher(pointPattern(:,1:2));
    [~,dmin] = knnsearch(kdtree,testPositions);
    F = mean(dmin(:) < distances, 1);
end