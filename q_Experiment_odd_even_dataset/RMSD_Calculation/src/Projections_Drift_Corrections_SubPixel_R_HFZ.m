function [projections1,projections2, ShiftM] = Projections_Drift_Corrections_SubPixel_R_HFZ(projections1,projections2, projections_cal, rangeY, rangeX)
fprintf('\nProjections Drift Corrections Rfacotr Version code: (multi-GPU) \n');

ShiftM = zeros(2,size(projections_cal,3));                                 % store the total shift
while true
    %% calculate the Factor
    FactorM = zeros(size(rangeY,2),size(rangeX,2),size(projections_cal,3));
    for tempy = 1:size(rangeY,2)
        for tempx = 1:size(rangeX,2)
            fprintf(['Position Y:',num2str(tempy),'/',num2str(size(rangeY,2)) ...
                     ' Position X:',num2str(tempx),'/',num2str(size(rangeX,2)),'\n']);
            projections_temp = zeros(size(projections1));
            parfor tempi = 1:size(projections_cal,3)
                projections_temp(:,:,tempi) = real(My_FourierShift(projections1(:,:,tempi),ShiftM(1,tempi)+rangeY(tempy),ShiftM(2,tempi)+rangeX(tempx)));
            end
            FactorM(tempy,tempx,:) = calculate_R_factor(projections_cal,projections_temp);
        end
    end
    clear projections_temp tempx tempy tempi
    
    %% find the positiion shift value
    flagShift = 0;
    for tempi = 1:size(projections_cal,3)
        FactorM_temp = FactorM(:,:,tempi);
        [~,linear_index] = min(FactorM_temp(:));
        % Convert the linear index to row and column indices
        [indY, indX] = ind2sub(size(FactorM_temp), linear_index);
        dy = rangeY(indY);ShiftM(1,tempi) = ShiftM(1,tempi)+dy;
        dx = rangeX(indX);ShiftM(2,tempi) = ShiftM(2,tempi)+dx;
        if abs(dy) == max(rangeY) || abs(dx) == max(rangeX)
            flagShift = 1;
        end
    end
    clear tempi FactorM_temp
    if flagShift == 0
        break
    end
end
for tempi = 1:size(projections_cal,3)
    % shift the projections
    projections1(:,:,tempi) = real(My_FourierShift(projections1(:,:,tempi),ShiftM(1,tempi),ShiftM(2,tempi)));
    projections2(:,:,tempi) = real(My_FourierShift(projections2(:,:,tempi),ShiftM(1,tempi),ShiftM(2,tempi)));
end
end