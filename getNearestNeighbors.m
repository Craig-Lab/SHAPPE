function nearest_neighbors = getNearestNeighbors(i,j,points,numNeigh)
    if(nargin < 4)
        numNeigh = 1;
    end
    points_i_index = find(points(:,end)==i);
    N_i = length(points_i_index);
    nearest_neighbors = zeros(numNeigh,N_i);
    for ni=1:N_i
        current_point_index = points_i_index(ni);
        current_x = points(current_point_index,1);
        current_y = points(current_point_index,2);
        points_j_index = find(points(:,end)==j);
        %Remove current point's index if i=j
        points_j_index = setdiff(points_j_index,current_point_index);
        cells_j_x = points(points_j_index,1);
        cells_j_y = points(points_j_index,2);
        all_distances = sqrt(((cells_j_x-current_x).^2)+((cells_j_y-current_y).^2));
        if(~isempty(points_j_index) && length(points_j_index)>=numNeigh)
            nearest_neighbors(1:numNeigh,ni) = mink(all_distances,numNeigh);
        elseif(~isempty(points_j_index))
            nearest_neighbors(1:length(points_j_index),ni) = sort(all_distances);
        end
    end
end