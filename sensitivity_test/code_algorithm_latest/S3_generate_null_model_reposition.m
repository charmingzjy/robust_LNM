clear; clc;

%% ===== Parameters =====
N_NULL = 10000;       % Number of null datasets to generate
MAX_TRIALS = 10000;   % Maximum placement attempts per lesion
rng('shuffle');

%% ===== Project paths =====

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
    sprintf('R2_null_model_random_indices_iter%d_reposition', N_NULL));

%% ===== Check inputs and create the output directory =====
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

% Remove empty rows and rows without a subject/lesion ID.
nonempty_rows = ~all(cellfun(@is_empty_cell, raw), 2);
data = raw(nonempty_rows, :);
data = data(~cellfun(@is_empty_cell, data(:, 1)), :);

if isempty(data)
    error('No lesion records were found in: %s', real_lesion_file);
end

% Remove a header row when the first NumInd value is not numeric.
if isnan(scalar_to_number(data{1, 2}))
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
valid_center_indices = find(brain_mask_linear);

if isempty(valid_center_indices)
    error('The brain mask contains no nonzero voxels.');
end

fprintf(['Brain mask loaded. Size = [%d %d %d], ' ...
    'nonzero voxels = %d.\n'], ...
    volume_size(1), volume_size(2), volume_size(3), ...
    numel(valid_center_indices));

%% ===== Preprocess the shape of each real lesion =====
real_lesions = struct( ...
    'Subid', {}, 'NumInd', {}, 'offsets', {}, 'orig_idx', {});

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

    % Store offsets from the first voxel so repositioning preserves the
    % original lesion shape and orientation.
    anchor = coordinates(1, :);
    offsets = coordinates - anchor;

    real_lesions(i).Subid = subject_ids{i};
    real_lesions(i).NumInd = num_indices(i);
    real_lesions(i).offsets = offsets;
    real_lesions(i).orig_idx = lesion_indices;
end

fprintf('All lesion shapes were preprocessed successfully.\n');

%% ===== Generate null datasets by random repositioning =====
for null_number = 1:N_NULL
    fprintf('\n=== Generating null %05d / %05d ===\n', ...
        null_number, N_NULL);

    null_data = struct('Subid', {}, 'NumInd', {}, 'lin_idx', {});

    for i = 1:n_lesions
        offsets = real_lesions(i).offsets;
        lesion_size = real_lesions(i).NumInd;

        [placed, new_indices] = place_lesion( ...
            offsets, lesion_size, valid_center_indices, ...
            volume_size, brain_mask_linear, MAX_TRIALS);

        if ~placed
            error(['Lesion %d (%s) could not be repositioned within ' ...
                '%d trials while generating null %05d.'], ...
                i, value_to_text(real_lesions(i).Subid, i), ...
                MAX_TRIALS, null_number);
        end

        null_data(i).Subid = real_lesions(i).Subid;
        null_data(i).NumInd = lesion_size;
        null_data(i).lin_idx = new_indices;
    end

    %% ===== Save the null dataset as Excel =====
    output_excel = fullfile( ...
        output_dir, sprintf('null_%05d.xlsx', null_number));

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

    %% ===== Save the null dataset as MAT =====
    output_mat = fullfile( ...
        output_dir, sprintf('null_%05d.mat', null_number));
    save(output_mat, 'null_data', '-v7.3');

    fprintf('Null %05d saved to:\n%s\n', null_number, output_mat);
end

fprintf('\nAll %d null datasets were generated successfully.\n', N_NULL);

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
