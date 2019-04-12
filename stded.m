function means = stded(gt, sim, w, h)
% Scaled Time Delay Embedding Distance
%
% Parameters
%
% gt                   (n,2) matrix of fixations from human scanpath
% sim                  (n,2) matrix of fixations from simulated scanpath
% w                    width of the stimulus image
% h                    height of the stimulus image
%
% Results:
%
% means                similarity measure
%
% Author               Vittorio Cuculo
% Date                 April 2019
%
% MATLAB porting from https://github.com/dariozanca/FixaTons/
%

% Normalize scanpaths coordinates
max_dim = max(w,h);
gt = gt / max_dim;
sim = sim / max_dim;

% Scanpath similarity is computed for all possible k
max_k = min(length(gt), length(sim));

similarities = zeros(1, max_k);
for k = 1:max_k
    s = tded(gt, sim, k, "Mean");
    similarities(k) = exp(-s);
    disp(similarities(k));
end

% Return the average similarity measure
means = mean(similarities);
end

function s = tded(gt, sim, k, mode)

% K must be smaller or equal then the length of scanpaths
if (length(gt) < k || length(sim) < k)
    error('ERROR: Too large value for the time-embedding vector dimension');
end

% Create time-embedding vectors for both scanpaths
n_gtvec = length(gt) - k + 1;
gt_vectors = cell(n_gtvec, 1);
for i = 1:n_gtvec
    gt_vectors{i} = gt(i:i+k-1,:);
end

n_simvec = length(sim) - k + 1;
sim_vectors = cell(n_simvec, 1);
for i = 1:n_simvec
    sim_vectors{i} = sim(i:i+k-1,:);
end

% Find the nearest vector between the k vectors from simulated scanpaths
% and the k vectors from human scanpaths
distances = zeros(n_simvec, 1);
for i = 1:n_simvec
    norms = zeros(n_gtvec, 1);
    for j = 1:n_gtvec
        d = norm(diag(pdist2(gt_vectors{j}, sim_vectors{i})));
        norms(j) = d;
    end
    distances(i) = min(norms) / k;
end

if mode == "Mean"
    s = mean(distances);
elseif mode == "Hausdorff"
    s = max(distances);
else
    error('ERROR: distance mode not defined.')
end
end