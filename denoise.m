function [vec_out] = denoise(vec_in)
%DENOISE Generates SAE (Surface of Active Events), do filtering (median
%filter), and outputs filtered event vectors.
vec_out = [];
SAE_temp = zeros(800, 1280);

for i = 1:length(vec_in)
    if (vec_in(i, 1) < 1 || vec_in(i, 2) < 1)
       continue 
    end
    SAE_temp(vec_in(i, 1), vec_in(i, 2)) = vec_in(i, 3);
end

SAE_size = size(SAE_temp);
SAE_temp_filtered = medfilt2(SAE_temp);

for i = 1:800
    for j = 1:1200
        if (SAE_temp_filtered(i, j) ~= 0)
            vec_temp = [i,j,SAE_temp_filtered(i, j)];
            vec_out = [vec_out; vec_temp];
        end
    end
end

end
