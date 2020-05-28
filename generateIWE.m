function [weightedIWE, evt_vec_warpped] = generateIWE(evt_vec, theta, P, ref_time)
%GENERATEIWE Computes an IWE (Image of Warpped Events) from a given event vector
%   Both weighted IWE and warpped event vector are outputs
size_theta = size(theta);
length_evt_vec = length(evt_vec);
cluster_num = size_theta(2);
weightedIWE = zeros(800,1280,cluster_num);
evt_vec_warpped = zeros(length_evt_vec,3,cluster_num);

for j = 1:cluster_num
    
    [row_warpped, col_warpped] = linearWarp(evt_vec(:,3) - ref_time, evt_vec(:,1), evt_vec(:,2), theta(1,j), theta(2,j));
    
    for k = 1:length_evt_vec % filter out out-of-bound events
        if (row_warpped(k) <= 1 || row_warpped(k) > 799 || ...
            col_warpped(k) <= 1 || col_warpped(k) > 1279)
            continue
        end
        
        temp_evt = [row_warpped(k), col_warpped(k), evt_vec(k,3)];
        evt_vec_warpped(k,:,j) = temp_evt; 
        
        temp_row = round(row_warpped(k)); % temp variables
        temp_col = round(col_warpped(k));
        
        weightedIWE(temp_row,temp_col,j) = weightedIWE(temp_row,temp_col,j) + P(k,j);       
    end
end

%% alternative computing method: compute warpped coordinates one by one, but slower
% for j = 1:length(theta)
%     for k = 1:length(evt_to_segment)
%         [row_warpped, col_warpped] = linearWarp(evt_to_segment(k,3) - ref_time, evt_to_segment(k,1), evt_to_segment(k,2), theta{j}(1), theta{j}(2));
%         evt_to_segment_warpped{j} = [evt_to_segment_warpped{j}; row_warpped, col_warpped, evt_to_segment(k,3)];
%         
%         if (row_warpped <= 1 || row_warpped > 799 || col_warpped <= 1 || col_warpped > 1279)
%             continue;
%         end
%         weightedIWE{j}(round(row_warpped), round(col_warpped)) = weightedIWE{j}(round(row_warpped), round(col_warpped)) + P(k,j);
%     end
% end
end

function [row_output, col_output] = linearWarp(delta_t, row_input, col_input, vel_row, vel_col)
%LINEARWARP Computes warpped event position, given event original coordinates, time interval, and linear velocities 
%   If delta_t, row_input and col_input are column vectors, then linearWarp also outputs
%   column vectors with the same dimensions.
row_output = row_input - vel_row.*delta_t;
col_output = col_input - vel_col.*delta_t;
end

