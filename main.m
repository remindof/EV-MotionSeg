%% Motion segmentation using motion compensation
% Jiang Rui@Celepixel   Apr-2020

clear all;

%% parameters - Offpixel
min_index = 268;
% noise test: 66, 182, 196, 234-236, 242-243, 256-257, 283-286, 301-303,
% 368-370, 384-387, 463-465, 517-518, 526-527
% one hand test: 20
% NOT SATISFACTORY: 234-236, 242-243, 283-286, 301-303, 463-465, 517-518, 526-527

max_index = 3000; % set the max index of csv files
num_clusters_in = 2; % set number of initial clusters
step_mu = 0.00000001; % set gradient ascent rate 0.00001 fine but slow
dist_thres = 1.0; % set dist threshold for delta function approximation, must be integer
thres_max_px_per_unit = 0.0015; % Off-pixel: set velocity threshold for initialization, cluster with larger velocity will be deemed as zero
thres_max_delta_t = 5000; % set the maximum delta_t considered in K-means clustering in function "initializeClusters"
thres_min_evt = 1000; % set minimum number of events in evt_vec for motion segmentation
thres_cost_non_assignment = 10; % set threshold for not assigning clusters to existing ones

%% main loop
cluster_dist_t = [];
variances = [];
variances_iter = [];
thetas_iter = [];

for i = min_index:max_index
    i
    %% data read-in and initialization
%     evt_vec = load(strcat(num2str(30000+i),'.csv'));
    evt_vec = load(strcat(num2str(i),'.csv'));
    evt_vec = denoise(evt_vec);
    if (length(evt_vec) < thres_min_evt)
        [cluster_dist_t] = [cluster_dist_t; 0];
        [variances] = [variances; 0 0];
        continue;
    end
    convergence = false;
    ref_time = evt_vec(floor(length(evt_vec)/2),3); % middle events used as the reference time
    
    %% Plotting
    figure(1)
    plot(-evt_vec(:,2), -evt_vec(:,1), '.');
    drawnow
    pause(0.1)
    %% initialize the cluster
    [theta, P, weightedIWE, evt_vec_warpped, cluster_center] = initializeClusters(evt_vec, num_clusters_in, ref_time, ...
        thres_max_px_per_unit, thres_max_delta_t);
    plotIWE(weightedIWE, i);
    %% use the cluster number to do motion segmentation
    iter_num = 1;
    while (convergence == false)
%         % recording data
        [variances_iter] = [variances_iter; imageVar(weightedIWE(:,:,1)) imageVar(weightedIWE(:,:,2))];
        [thetas_iter] = [thetas_iter; theta(1,1) theta(2,1) theta(1,2) theta(2,2)];
        
%         [variances_iter] = [variances_iter; imageVar(weightedIWE(:,:,1)) ];
%         [thetas_iter] = [thetas_iter; theta(1,1) theta(2,1) ];
        
        % updating params
        [P] = updateAssignments(evt_vec_warpped, weightedIWE, P);
        [weightedIWE, evt_vec_warpped] = generateIWE(evt_vec, theta, P, ref_time);
        %[theta, weightedIWE, evt_vec_warpped] = updateMotionParam(evt_vec, weightedIWE, evt_vec_warpped, theta, P, ref_time, step_mu, dist_thres);
        
        iter_num = iter_num + 1;
        
        plotIWE(weightedIWE, i);
        if (iter_num > 50)
            convergence = true;
        end
    end
    %% plotting
    plotIWE(weightedIWE, i);
%     plotProbMap(evt_vec, P)

    %% for testing only
    [cluster_dist_t] = [cluster_dist_t; norm(cluster_center(:,1)-cluster_center(:,2))];
    [variances] = [variances; imageVar(weightedIWE(:,:,1)) imageVar(weightedIWE(:,:,2))];
        
end

%% plotting functions
function plotRawEvts(evt_vec)
figure(11)
plot(evt_vec(:,2), -evt_vec(:,1), '.', 'MarkerSize',1);
xlim([0 1280])
ylim([-800 0])
axis equal
drawnow
pause(0.01)
end

function plotProbMap(evt_vec, P)
size_P = size(P); % row is the event number, col is the cluster number
for i = 1:size_P(2)
    figure(30+i);
    scatter(evt_vec(:,2),-evt_vec(:,1), 1, P(:,i));
    title(strcat('Cluster ',num2str(i)));
    colorbar;colormap('winter');caxis([0,1]);
end
end

function plotIWE(weightedIWE, index)
size_weightedIWE = size(weightedIWE);
if (length(size_weightedIWE) == 3)
for i = 1:size_weightedIWE(3)
    figure(90+i);
    imshow(weightedIWE(:,:,i));
    saveas(gcf, strcat(num2str(index,'%03d'), num2str(i), '.png'));
end
else
    figure(90);
    imshow(weightedIWE(:,:));
end
end







