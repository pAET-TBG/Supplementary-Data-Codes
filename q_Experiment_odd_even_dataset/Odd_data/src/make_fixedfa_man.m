function fixedfa =  make_fixedfa_man(sizeX, Res, Z_arr)
    finalvol_summed = zeros(sizeX);                                        % initialize the finalvol_summed (the projections)
        
    kx = 1:size(finalvol_summed,1);                                        % the vector for column in reciprocal space
    ky = 1:size(finalvol_summed,2);                                        % the vector for row in reciprocal space

    MultF_X = 1/(length(kx)*Res);                                          % the pixel size for column in reciprocal space
    MultF_Y = 1/(length(ky)*Res);                                          % the pixel size for row in reciprocal space

    CentPos = round((size(finalvol_summed)+1)/2);                          % the center of column CentPos(1) and row CentPos(2)
    [KX, KY] = ndgrid((kx-CentPos(1))*MultF_X,(ky-CentPos(2))*MultF_Y);    % the grid for the reciprocal space KX is for column, KY is for row   
    q2 = KX.^2 + KY.^2 ;                                                   % the square of distance between each positions in the grid with center
    clear KX KY KZ

    fixedfa_arr = zeros(numel(Z_arr), numel(q2));                          % initialize the fixedfa_arr marix of (3, sizeX(1)*sizeX(2))
    for i = 1:numel(Z_arr)
        fixedfa_arr(i,:) = fatom_vector(sqrt(q2),Z_arr(i));                % calculate each value in fixedfa_arr
    end
    fixedfa = mean(fixedfa_arr,1);                                         % average the column value matrix from (3, sizeX(1)*sizeX(2)) to (1,  sizeX(1)*sizeX(2))
end