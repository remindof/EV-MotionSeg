function [theta_out, weightedIWE_out, evt_vec_warpped_out] = updateMotionParam(evt_vec, ...
    weightedIWE_in, evt_vec_warpped_in, theta_in, P, ref_time, step_mu, dist_thres)
%UPDATEMOTIONPARAM Update motion parameter (linear optical flow/velocity in
%this case) towards maximum the sum of contrasts of all clusters using gradient ascent
size_theta_in = size(theta_in);
cluster_num = size_theta_in(2);

theta_out = zeros(2,cluster_num);

for j = 1:cluster_num
    grad_theta_j = findGrad(evt_vec_warpped_in(:,:,j), weightedIWE_in(:,:,j), P(:,j), dist_thres, ref_time);
    theta_out(:,j) = theta_in(:,j) + step_mu .* grad_theta_j;
end

% generate new IWE using raw event vector, probability matrix, and updated motion parameters
[weightedIWE_out, evt_vec_warpped_out] = generateIWE(evt_vec, theta_out, P, ref_time);

end

function [grad_theta_j] = findGrad(evt_vec_warpped_j, weightedIWE_j, P_j, dist_thres, ref_time)
row_min = round(max( [min(evt_vec_warpped_j(:,1)),1] ));
row_max = round(min( [max(evt_vec_warpped_j(:,1)),799] ));
col_min = round(max( [min(evt_vec_warpped_j(:,2)),1] ));
col_max = round(min( [max(evt_vec_warpped_j(:,2)),1279] ));
grad_I_j_row = zeros(800,1280);
grad_I_j_col = zeros(800,1280);
sum_grad_I_j_row = 0;
sum_grad_I_j_col = 0;
mu_j = imageMean(weightedIWE_j(row_min:row_max, col_min:col_max));
%delta_t = evt_vec_warpped_j(:,3) - ref_time;

% optimized code: no performance improved
%         dist_row = row - evt_vec_warpped_j(:,1);
%         dist_col = col - evt_vec_warpped_j(:,2); % x - x'(kj)
%         grad_I_j_row(row,col) = sum(P_j(:) .* findGradDelta(dist_row, dist_thres) .* (-delta_t)); % eqn(1)
%         grad_I_j_col(row,col) = sum(P_j(:) .* findGradDelta(dist_col, dist_thres) .* (-delta_t));
% optimized code: end

for row = row_min:row_max
    for col = col_min:col_max
        for k = 1:length(evt_vec_warpped_j)
            delta_t_k = evt_vec_warpped_j(k,3) - ref_time;
            dist_row = row - evt_vec_warpped_j(k,1);
            dist_col = col - evt_vec_warpped_j(k,2); % x - x'(kj)
            if ( (abs(dist_row)>dist_thres) ) % || (dist_row==0) ) % for row
            
            else
                grad_I_j_row(row,col) = grad_I_j_row(row,col) + P_j(k) * findGradDelta(dist_row) * (delta_t_k); % eqn(1)
            end
            if ( (abs(dist_col)>dist_thres) ) % || (dist_col==0) ) % for col
            
            else
                grad_I_j_col(row,col) = grad_I_j_col(row,col) + P_j(k) * findGradDelta(dist_col) * (delta_t_k);
            end
        end
        sum_grad_I_j_row = sum_grad_I_j_row + grad_I_j_row(row,col);
        sum_grad_I_j_col = sum_grad_I_j_col + grad_I_j_col(row,col);
    end 
end
avg_grad_I_j_row = sum_grad_I_j_row./((row_max - row_min)*(col_max - col_min)); % eqn(2)
avg_grad_I_j_col = sum_grad_I_j_col./((row_max - row_min)*(col_max - col_min));

sum_x_row = 0;
sum_x_col = 0;
for row = row_min:row_max
    for col = col_min:col_max
        sum_x_row = sum_x_row + (weightedIWE_j(row,col) - mu_j) * (grad_I_j_row(row,col) - avg_grad_I_j_row);
        sum_x_col = sum_x_col + (weightedIWE_j(row,col) - mu_j) * (grad_I_j_col(row,col) - avg_grad_I_j_col);
    end
end
grad_theta_j = (2/((row_max - row_min)*(col_max - col_min))).*[sum_x_row; sum_x_col];
end

% function [grad_delta] = findGradDelta(dist, dist_thres)
% grad_delta = -0.1 .* (abs(dist(:))<=dist_thres) .* dist .^ 3;
% end

function [grad_delta] = findGradDelta(dist)
%FINDGRADDELTA Returns gradient of approximated Dirac delta function
grad_delta = -dist^3;
end

function [img_mean] = imageMean(img)
%IMAGEMEAN Computes the mean of an image
img_mean = mean(img(:));
end