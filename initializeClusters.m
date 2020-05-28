function [theta, P, weightedIWE, evt_vec_warpped, cluster_center] ...
    = initializeClusters(evt_vec, num_clusters_in, ref_time, thres_max_px_per_unit, thres_max_delta_t)
%INITIALIZECLUSTERS Create clusters, get initial parameters done
P = (1./num_clusters_in).*ones(length(evt_vec), num_clusters_in); % dimension = evt num * cluster num
theta = zeros(2,num_clusters_in); % each col vector in theta has the dimension of 2 (row vel and col vel)
cluster_center = zeros(2,num_clusters_in); % same format as theta

% automatically initialize theta
[theta, cluster_center] = kMeansClusterOpticalFlow(evt_vec, num_clusters_in, ...
    thres_max_px_per_unit, thres_max_delta_t);

[weightedIWE, evt_vec_warpped] = generateIWE(evt_vec, theta, P, ref_time);
end

function [theta, cluster_center] = kMeansClusterOpticalFlow(evt_vec, num_clusters_in, ...
    thres_max_px_per_unit, thres_max_delta_t)

theta = zeros(2,num_clusters_in);
optical_flow_vec = [];
SAE = zeros(800,1280);

for i = 1:length(evt_vec) % generating SAE, find optical flow
    if (SAE(evt_vec(i,1), evt_vec(i,2)) < evt_vec(i,3))
        SAE(evt_vec(i,1), evt_vec(i,2)) = evt_vec(i,3);
    end
end

diff_row_SAE = diff(SAE,1,1);
diff_col_SAE = diff(SAE,1,2);

for row = 1:800-1
    for col = 1:1280-1
        if (SAE(row,col)==0 || diff_row_SAE(row,col)==0 || diff_col_SAE(row,col)==0 || diff_row_SAE(row,col)==diff_col_SAE(row,col) ...
         || SAE(row+1,col)==0 || SAE(row,col+1)==0 || abs(diff_row_SAE(row,col))>thres_max_delta_t || abs(diff_col_SAE(row,col))>thres_max_delta_t )
            continue
        end
        optical_flow_vec = [optical_flow_vec; row col diff_row_SAE(row,col) diff_col_SAE(row,col)];
    end
end
[idx,cluster_center] = kmeans([optical_flow_vec(:,3), optical_flow_vec(:,4)], num_clusters_in, 'Replicates',2);
%[idx,cluster_center] = kmeans([optical_flow_vec(:,3), optical_flow_vec(:,4)], num_clusters_in, 'Distance','cityblock','Replicates',2);
cluster_center = cluster_center';
theta = 1./cluster_center;

for i = 1:num_clusters_in
    if (theta(1,i) == inf || abs(theta(1,i)) > thres_max_px_per_unit)
        theta(1,i) = 0;
    end
    if (theta(2,i) == inf || abs(theta(2,i)) > thres_max_px_per_unit)
        theta(2,i) = 0;
    end
end

%% plot K-means results
figure(200);
plot(optical_flow_vec(idx==1,3),optical_flow_vec(idx==1,4),'r.','MarkerSize',5)

hold on
plot(optical_flow_vec(idx==2,3),optical_flow_vec(idx==2,4),'b.','MarkerSize',5)
plot(cluster_center(1,:),cluster_center(2,:),'kx','MarkerSize',15,'LineWidth',3)
legend('Cluster 1','Cluster 2','Centroids','Location','NW')
title 'Cluster Assignments and Centroids'
hold off

end