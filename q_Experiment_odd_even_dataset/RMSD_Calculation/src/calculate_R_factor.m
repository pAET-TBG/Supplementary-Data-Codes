function R_factor = calculate_R_factor(projections_cal,projections_exp)
R_factor = zeros(size(projections_cal,3),1);
for tempi = 1:size(projections_cal,3)
    R_factor(tempi) = sum(sum(abs(projections_cal(:,:,tempi)-projections_exp(:,:,tempi))))/sum(sum(abs(projections_exp(:,:,tempi))));
end
end

