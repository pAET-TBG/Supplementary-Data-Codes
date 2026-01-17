function twist_angle = calculate_twist_angle(atompos, theta0, searchRange, atomDist, flag)
% atompos:     the atom coordinates of the bilayer graphene
% theta0:      roughpy initial guess of the twisting angle
% searchRange: the searching range of the twisting angle
% flag:        whether to visualize the matching result after twisting back

atomposup = atompos(1:2,atompos(3,:)<0);
atomposdown = atompos(1:2,atompos(3,:)>0);
% rotation matrix
R0 = [cosd(theta0) -sind(theta0); sind(theta0) cosd(theta0)];
atomposup = R0*atomposup;

thetaVector = searchRange;
dirVector = zeros(1,size(thetaVector,2));
count = 1;
%% match the upper and lower atom 
Mag_shift = 1;
shift=[0 0]';   % the shift used to do the subpixel shift to align the atom position
while Mag_shift>1e-5 % if the total difference is larger than 1e-6 (angstrom) then continue to do the alignment
    atomposup = atomposup  - shift;
    difarr = [];    % array for store the difference between two tracing result
    disarr = [];    % array for store the distance between two tracing result
    count_arr1 = [];% store the index for the first atom result, which is match with the second atom tracing result
    count_arr2 = [];% store the index for the second atom result, which is match with the first atom tracing result
    for i=1:size(atomposdown,2)
        dif=(atomposup-atomposdown(:,i));         % calculate all difference from the first result with i'th position in the second result (angstrom)
        dis=sqrt(sum(dif.^2,1));                  % calculate the distance (angstrom)
        [dis,ind]=min(dis);                       % obatin the minimum distance and corresponding index
        if dis <= atomDist                        % if the minimum distance is smaller than the threshold then store the information
            difarr=[difarr dif(:,ind)];
            disarr = [disarr dis];
        end
    end
    shift=mean(difarr,2);                  % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
    Mag_shift=sum(abs(mean(difarr,2)));          % calculate the total difference
end

for tempTheta = thetaVector
    tic
    % rotation matrix
    R0 = [cosd(tempTheta) -sind(tempTheta); sind(tempTheta) cosd(tempTheta)];
    tempatomposup = R0*atomposup;
    Mag_shift = 1;
    shift=[0 0]';   % the shift used to do the subpixel shift to align the atom position
    while Mag_shift>1e-5 % if the total difference is larger than 1e-6 (angstrom) then continue to do the alignment
        tempatomposup = tempatomposup  - shift;
        difarr = [];    % array for store the difference between two tracing result
        disarr = [];    % array for store the distance between two tracing result
        count_arr1 = [];% store the index for the first atom result, which is match with the second atom tracing result
        count_arr2 = [];% store the index for the second atom result, which is match with the first atom tracing result
        for i=1:size(atomposdown,2)
            dif=(tempatomposup-atomposdown(:,i));         % calculate all difference from the first result with i'th position in the second result (angstrom)
            dis=sqrt(sum(dif.^2,1));                  % calculate the distance (angstrom)
            [dis,ind]=min(dis);                       % obatin the minimum distance and corresponding index
            if dis <= atomDist                        % if the minimum distance is smaller than the threshold then store the information
                difarr=[difarr dif(:,ind)];
                disarr = [disarr dis];
            end
        end
        shift=mean(difarr,2);                  % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
        Mag_shift=sum(abs(mean(difarr,2)));          % calculate the total difference
    end
    dirVector(count) = mean(disarr);
    count = count + 1;
    toc
end
[~, index] = min(dirVector);
theta1 = thetaVector(index);
R0 = [cosd(theta1) -sind(theta1); sind(theta1) cosd(theta1)];
atomposup = R0*atomposup;

Mag_shift = 1;
meandisarr = inf;
shift=[0 0]';   % the shift used to do the subpixel shift to align the atom position
while Mag_shift>1e-5 % if the total difference is larger than 1e-6 (angstrom) then continue to do the alignment
    atomposup = atomposup  - shift;
    difarr = [];    % array for store the difference between two tracing result
    disarr = [];    % array for store the distance between two tracing result
    count_arr1 = [];% store the index for the first atom result, which is match with the second atom tracing result
    count_arr2 = [];% store the index for the second atom result, which is match with the first atom tracing result
    for i=1:size(atomposdown,2)
        dif=(atomposup-atomposdown(:,i));         % calculate all difference from the first result with i'th position in the second result (angstrom)
        dis=sqrt(sum(dif.^2,1));                  % calculate the distance (angstrom)
        [dis,ind]=min(dis);                       % obatin the minimum distance and corresponding index
        if dis <= atomDist                        % if the minimum distance is smaller than the threshold then store the information
            difarr=[difarr dif(:,ind)];
            disarr = [disarr dis];
        end
    end
    shift=mean(difarr,2);                  % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
    Mag_shift=sum(abs(mean(difarr,2)));          % calculate the total difference
end

if flag == 1
    % visualize the match region if the we correctly calculate the twisting
    % angle the upper layer and the lower layer will match perfectly after
    % twist the upper layer back
    figure();plot(atomposdown(1,:),atomposdown(2,:),'b.');hold on;
    plot(atomposup(1,:),atomposup(2,:),'r.');hold off;
end

twist_angle = theta1+theta0;

end