function [P_out] = updateAssignments(evt_vec_warpped, weightedIWE, P_in)
%UPDATEASSIGNMENTS Updates probabilities that an event belongs to a cluster
P_out = P_in;
size_evt_vec_warpped = size(evt_vec_warpped);
cluster_num = length(P_in(1,:));

for k = 1:size_evt_vec_warpped(1) % for each event in the first cluster
    sum_c = 0;
    P_row_temp = zeros(1,cluster_num); % for each cluster
    
    for j = 1:cluster_num
        
        row_warpped = round(evt_vec_warpped(k,1,j));
        col_warpped = round(evt_vec_warpped(k,2,j));
        
        if (row_warpped == 0 || col_warpped == 0)
            continue;
        end
        
        if (weightedIWE(row_warpped,col_warpped,j) ~= 0) % update probability
            sum_c = sum_c + weightedIWE(row_warpped,col_warpped,j);
            P_row_temp(1,j) = weightedIWE(row_warpped,col_warpped,j);
        end
    end
    
    if (sum_c ~= 0)
        P_out(k,:) = P_row_temp ./ sum_c;
    end
end
end

