clear; clc;

%% ===== Parameters =====
N_NULL = 10000;       % Number of qualified null datasets to retain
MAX_TRIALS = 1000;    % Maximum placement attempts per lesion
N_OPT_ITER = 1000;    % Simulated annealing steps per candidate dataset

% Reduce T0 if almost every proposed move is accepted. Increase it slightly
% if the optimizer frequently becomes trapped.
T0 = 1e-4;            % Initial simulated annealing temperature
ALPHA = 0.997;        % Cooling coefficient

P_THRESHOLD = 0.05;   % Accept only if the KS-test p-value is greater than this
rng('shuffle');

%% ===== Project paths =====
% Expected project structure:
% null_model_pipeline/
%   code_pipeline/
%     S3_generate_null_model_annealing.m
%   disease_data/
%     Amnesia_mask_voxel_indices_v2.xlsx
%   template/
%     norm_8mm.nii.gz
%   output/
%     R2_null_model_random_indices_iter10000_sortedfreq/
%
% Build every path relative to the location of this script.
script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(script_dir);

disease_data_dir = fullfile(project_dir, 'disease_data');
template_dir = fullfile(project_dir, 'template');
project_output_dir = fullfile(project_dir, 'output');

real_lesion_file = fullfile( ...
    disease_data_dir, 'Amnesia_mask_voxel_indices_v2.xlsx');
brain_mask_file = fullfile(template_dir, 'norm_8mm.nii.gz');

output_dir = fullfile( ...
    project_output_dir, ...
    sprintf('R2_null_model_random_indices_iter%d_sortedfreq', N_NULL));

%% ===== Check input files and create the output directory =====
if ~isfile(real_lesion_file)
    error('Real lesion index file not found: %s', real_lesion_file);
end

if ~isfile(brain_mask_file)
    error('Brain mask file not found: %s', brain_mask_file);
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% ===== Load real lesion information =====
% Expected Excel format:
% Column 1: subject/lesion ID
% Column 2: number of lesion voxels
% Column 3 onward: voxel linear indices
raw = readcell(real_lesion_file);

if size(raw, 2) < 3
    error('The lesion index file must contain at least three columns.');
end

% Remove completely empty rows and rows without a subject/lesion ID.
nonempty_rows = ~all(cellfun(@is_empty_cell, raw), 2);
data = raw(nonempty_rows, :);
data = data(~cellfun(@is_empty_cell, data(:, 1)), :);

if isempty(data)
    error('No lesion records were found in: %s', real_lesion_file);
end

% Automatically remove a header row when column 2 is not numeric.
first_num_ind = scalar_to_number(data{1, 2});
if isnan(first_num_ind)
    data(1, :) = [];
end

if isempty(data)
    error('No lesion records remained after removing the header row.');
end

n_lesions = size(data, 1);
subject_ids = data(:, 1);
num_indices = zeros(n_lesions, 1);

for i = 1:n_lesions
    num_indices(i) = scalar_to_number(data{i, 2});

    if isnan(num_indices(i)) || num_indices(i) < 1 || ...
            num_indices(i) ~= floor(num_indices(i))
        error('Invalid NumInd value in lesion row %d.', i);
    end
end

fprintf('Loaded %d lesions from:\n%s\n', n_lesions, real_lesion_file);

%% ===== Load the brain mask =====
if exist('MRIread', 'file') ~= 2
    error(['MRIread was not found on the MATLAB path. Add FreeSurfer ' ...
        'MATLAB utilities before running this script.']);
end

mri = MRIread(brain_mask_file);
brain_mask = mri.vol ~= 0;
volume_size = size(brain_mask);

if numel(volume_size) ~= 3
    error('The brain mask must be a three-dimensional volume.');
end

num_volume_voxels = numel(brain_mask);
brain_mask_linear = brain_mask(:);

% The same brain-mask voxels are used both as candidate anchor locations
% and as the region for frequency-distribution matching.
valid_indices = find(brain_mask_linear);
valid_center_indices = valid_indices;

if isempty(valid_indices)
    error('The brain mask contains no nonzero voxels.');
end

fprintf(['Brain mask loaded. Size = [%d %d %d], ' ...
    'nonzero voxels = %d.\n'], ...
    volume_size(1), volume_size(2), volume_size(3), ...
    numel(valid_indices));

%% ===== Preprocess lesion shapes and the true frequency map =====
real_lesions = struct( ...
    'Subid', {}, 'NumInd', {}, 'offsets', {}, 'orig_idx', {});
true_frequency = zeros(num_volume_voxels, 1);

for i = 1:n_lesions
    subject_id_text = value_to_text(subject_ids{i}, i);
    expected_count = num_indices(i);

    lesion_cells = data(i, 3:end);

    lesion_indices = [];
    for k = 1:numel(lesion_cells)
        voxel_index = scalar_to_number(lesion_cells{k});

        if ~isnan(voxel_index)
            lesion_indices(end + 1, 1) = voxel_index; %#ok<SAGROW>
        end
    end

    lesion_indices = unique(lesion_indices, 'stable');

    if isempty(lesion_indices)
        error('No voxel indices were found for lesion %d (%s).', ...
            i, subject_id_text);
    end

    if numel(lesion_indices) ~= expected_count
        warning(['Lesion %d (%s): NumInd = %d, but %d unique indices ' ...
            'were read. NumInd was updated.'], ...
            i, subject_id_text, expected_count, numel(lesion_indices));
        num_indices(i) = numel(lesion_indices);
    end

    if any(lesion_indices < 1 | lesion_indices > num_volume_voxels | ...
            lesion_indices ~= floor(lesion_indices))
        error('Lesion %d (%s) contains invalid voxel indices.', ...
            i, subject_id_text);
    end

    if ~all(brain_mask_linear(lesion_indices))
        warning('Lesion %d (%s) contains voxels outside the brain mask.', ...
            i, subject_id_text);
    end

    [x, y, z] = ind2sub(volume_size, lesion_indices);
    coordinates = [x(:), y(:), z(:)];

    % Preserve the lesion shape through integer translation.
    anchor = coordinates(1, :);
    offsets = coordinates - anchor;

    real_lesions(i).Subid = subject_ids{i};
    real_lesions(i).NumInd = num_indices(i);
    real_lesions(i).offsets = offsets;
    real_lesions(i).orig_idx = lesion_indices;

    true_frequency(lesion_indices) = ...
        true_frequency(lesion_indices) + 1;
end

true_ratio = true_frequency / n_lesions;
true_values_valid = true_ratio(valid_indices);
true_sorted_valid = sort(true_values_valid, 'descend');

fprintf('The true lesion-frequency distribution has been computed.\n');

%% ===== Generate null datasets =====
num_accepted = 0;
attempt = 0;

while num_accepted < N_NULL
    attempt = attempt + 1;

    fprintf('\n=== Attempt %05d | Accepted %05d / %05d ===\n', ...
        attempt, num_accepted, N_NULL);

    %% ----- 1. Generate an initial valid null dataset -----
    null_data = struct('Subid', {}, 'NumInd', {}, 'lin_idx', {});
    current_frequency = zeros(num_volume_voxels, 1);
    initial_solution_valid = true;

    for i = 1:n_lesions
        offsets = real_lesions(i).offsets;
        lesion_size = real_lesions(i).NumInd;

        [placed, new_indices] = place_lesion( ...
            offsets, lesion_size, valid_center_indices, volume_size, ...
            brain_mask_linear, MAX_TRIALS);

        if ~placed
            warning(['Initial placement failed for lesion %d (%s). ' ...
                'The candidate null dataset was discarded.'], ...
                i, value_to_text(real_lesions(i).Subid, i));
            initial_solution_valid = false;
            break;
        end

        null_data(i).Subid = real_lesions(i).Subid;
        null_data(i).NumInd = lesion_size;
        null_data(i).lin_idx = new_indices;

        current_frequency(new_indices) = ...
            current_frequency(new_indices) + 1;
    end

    if ~initial_solution_valid
        continue;
    end

    current_ratio = current_frequency / n_lesions;
    current_loss = sorted_frequency_loss( ...
        current_ratio, true_sorted_valid, valid_indices);

    %% ----- 2. Optimize the frequency distribution -----
    temperature = T0;
    num_move_attempts = 0;
    num_moves_accepted = 0;

    for iteration = 1:N_OPT_ITER
        lesion_number = randi(n_lesions);

        old_indices = null_data(lesion_number).lin_idx;
        offsets = real_lesions(lesion_number).offsets;
        lesion_size = real_lesions(lesion_number).NumInd;

        [placed, new_indices] = place_lesion( ...
            offsets, lesion_size, valid_center_indices, volume_size, ...
            brain_mask_linear, MAX_TRIALS);

        if ~placed
            temperature = temperature * ALPHA;
            continue;
        end

        num_move_attempts = num_move_attempts + 1;

        proposed_frequency = current_frequency;
        proposed_frequency(old_indices) = ...
            proposed_frequency(old_indices) - 1;
        proposed_frequency(new_indices) = ...
            proposed_frequency(new_indices) + 1;

        proposed_ratio = proposed_frequency / n_lesions;
        proposed_loss = sorted_frequency_loss( ...
            proposed_ratio, true_sorted_valid, valid_indices);

        if proposed_loss < current_loss
            accept_move = true;
        else
            acceptance_probability = exp( ...
                -(proposed_loss - current_loss) / ...
                max(temperature, eps));
            accept_move = rand < acceptance_probability;
        end

        if accept_move
            current_frequency = proposed_frequency;
            current_ratio = proposed_ratio;
            current_loss = proposed_loss;
            null_data(lesion_number).lin_idx = new_indices;
            num_moves_accepted = num_moves_accepted + 1;
        end

        temperature = temperature * ALPHA;
    end

    %% ----- 3. Compare the true and null frequency distributions -----
    null_values_valid = current_ratio(valid_indices);
    [ks_reject, ks_p_value] = kstest2( ...
        true_values_valid, null_values_valid, ...
        'Alpha', P_THRESHOLD);

    if num_move_attempts > 0
        move_acceptance_rate = ...
            num_moves_accepted / num_move_attempts;
    else
        move_acceptance_rate = NaN;
    end

    fprintf(['Final sorted-frequency RMSE = %.8f | KS p = %.6g | ' ...
        'h = %d | move acceptance rate = %.3f\n'], ...
        current_loss, ks_p_value, ks_reject, move_acceptance_rate);

    if isnan(ks_p_value) || ks_reject == 1 || ...
            ks_p_value <= P_THRESHOLD
        fprintf('Rejected: the KS test was significant or invalid.\n');
        continue;
    end

    %% ----- 4. Accept and save the null dataset -----
    num_accepted = num_accepted + 1;
    fprintf('Accepted as null_%05d.\n', num_accepted);

    output_excel = fullfile( ...
        output_dir, sprintf('null_%05d.xlsx', num_accepted));

    max_lesion_size = max(num_indices);
    output_table = cell(n_lesions + 1, 2 + max_lesion_size);
    output_table(1, 1:3) = {'Subid', 'NumInd', 'Ind'};

    for i = 1:n_lesions
        lesion_size = null_data(i).NumInd;
        lesion_indices = null_data(i).lin_idx(:);

        output_table{i + 1, 1} = null_data(i).Subid;
        output_table{i + 1, 2} = lesion_size;
        output_table(i + 1, 3:(2 + lesion_size)) = ...
            num2cell(lesion_indices');
    end

    writecell(output_table, output_excel);

    output_mat = fullfile( ...
        output_dir, sprintf('null_%05d.mat', num_accepted));

    null_ratio = current_ratio;
    final_loss = current_loss;
    final_acceptance_rate = move_acceptance_rate;
    p_ks = ks_p_value;
    h_ks = ks_reject;

    save(output_mat, ...
        'null_data', 'null_ratio', 'p_ks', 'h_ks', ...
        'valid_indices', 'final_loss', 'final_acceptance_rate', ...
        'true_sorted_valid', '-v7.3');
end

fprintf('\nAll %d qualified null datasets were generated.\n', N_NULL);

%% ===== Local functions =====
function [placed, linear_indices] = place_lesion( ...
        offsets, lesion_size, valid_centers, volume_size, ...
        brain_mask_linear, max_trials)

    placed = false;
    linear_indices = [];

    for trial = 1:max_trials
        center_linear = valid_centers(randi(numel(valid_centers)));
        [center_x, center_y, center_z] = ...
            ind2sub(volume_size, center_linear);

        new_coordinates = offsets + repmat( ...
            [center_x, center_y, center_z], lesion_size, 1);

        x = new_coordinates(:, 1);
        y = new_coordinates(:, 2);
        z = new_coordinates(:, 3);

        if any(x < 1 | x > volume_size(1) | ...
                y < 1 | y > volume_size(2) | ...
                z < 1 | z > volume_size(3))
            continue;
        end

        candidate_indices = sub2ind(volume_size, x, y, z);

        if ~all(brain_mask_linear(candidate_indices))
            continue;
        end

        if numel(unique(candidate_indices)) ~= lesion_size
            continue;
        end

        linear_indices = candidate_indices(:);
        placed = true;
        return;
    end
end

function loss = sorted_frequency_loss( ...
        null_ratio, true_sorted_valid, valid_indices)

    null_sorted_valid = sort(null_ratio(valid_indices), 'descend');
    differences = null_sorted_valid - true_sorted_valid;
    loss = sqrt(mean(differences .^ 2));
end

function number = scalar_to_number(value)
    if isnumeric(value) && isscalar(value)
        number = double(value);
    elseif islogical(value) && isscalar(value)
        number = double(value);
    elseif ischar(value) || (isstring(value) && isscalar(value))
        number = str2double(value);
    else
        number = NaN;
    end
end

function is_empty = is_empty_cell(value)
    if isempty(value)
        is_empty = true;
    elseif ismissing(value)
        is_empty = true;
    elseif ischar(value)
        is_empty = isempty(strtrim(value));
    elseif isstring(value) && isscalar(value)
        is_empty = strlength(strtrim(value)) == 0;
    elseif isnumeric(value) && isscalar(value)
        is_empty = isnan(value);
    else
        is_empty = false;
    end
end

function text = value_to_text(value, fallback_number)
    if iscell(value) && isscalar(value)
        value = value{1};
    end

    if isnumeric(value) && isscalar(value)
        text = num2str(value);
    elseif ischar(value)
        text = value;
    elseif isstring(value) && isscalar(value)
        text = char(value);
    elseif iscategorical(value) && isscalar(value)
        text = char(value);
    else
        text = sprintf('unknownID_%04d', fallback_number);
    end
end
