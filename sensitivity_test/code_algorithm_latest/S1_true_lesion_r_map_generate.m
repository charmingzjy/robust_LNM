clear; clc;

%% ===== Project paths =====
script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(script_dir);

disease_data_dir = fullfile(project_dir, 'disease_data');
template_dir = fullfile(project_dir, 'template');
project_output_dir = fullfile(project_dir, 'output');

lesion_indices_file = fullfile( ...
    disease_data_dir, 'Amnesia_mask_voxel_indices_v2.xlsx');
valid_indices_file = fullfile( ...
    template_dir, 'valid_indices_standard_wh_MNI8mm_32_32_32.mat');
mean_fc_file = fullfile(template_dir, 'rho_z_mean.mat');

output_dir = fullfile( ...
    project_output_dir, 'R0_true_lesion_r_maps_from_meanFC');

%% ===== Check input files =====
required_files = {
    lesion_indices_file
    valid_indices_file
    mean_fc_file
};

for i_file = 1:numel(required_files)
    if ~isfile(required_files{i_file})
        error('Required input file not found: %s', required_files{i_file});
    end
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% ===== Load lesion voxel indices =====
% Expected Excel format:
% Column 1: subject/lesion ID
% Column 2: number of lesion voxels
% Column 3 onward: voxel linear indices
data = readtable(lesion_indices_file, 'Range', 'A1');

fprintf('Loaded %d lesions from:\n%s\n', height(data), lesion_indices_file);

%% ===== Load valid voxel indices =====
valid_indices_data = load(valid_indices_file);

if ~isfield(valid_indices_data, 'valid_indices')
    error('Variable "valid_indices" was not found in: %s', ...
        valid_indices_file);
end

valid_indices = valid_indices_data.valid_indices(:);

if isempty(valid_indices)
    error('Variable "valid_indices" is empty.');
end

fprintf('Loaded %d valid voxel indices.\n', numel(valid_indices));

%% ===== Load the GSP1000 mean functional-connectivity matrix =====
mean_fc_data = load(mean_fc_file);

if isfield(mean_fc_data, 'rho_z')
    rho_z_mean = mean_fc_data.rho_z;
elseif isfield(mean_fc_data, 'rho_z_mean')
    rho_z_mean = mean_fc_data.rho_z_mean;
else
    error(['Neither "rho_z" nor "rho_z_mean" was found in: %s'], ...
        mean_fc_file);
end

if size(rho_z_mean, 1) ~= numel(valid_indices)
    error(['The number of rows in the mean FC matrix (%d) does not ' ...
        'match the number of valid indices (%d).'], ...
        size(rho_z_mean, 1), numel(valid_indices));
end

fprintf('Mean FC matrix loaded. Size = %d x %d.\n', ...
    size(rho_z_mean, 1), size(rho_z_mean, 2));

%% ===== Generate one mean FC map for each lesion =====
for i = 1:height(data)

    tic;

    % Convert the subject/lesion ID to text.
    sub_id_raw = data{i, 1};

    if iscell(sub_id_raw)
        sub_id_raw = sub_id_raw{1};
    end

    if isnumeric(sub_id_raw)
        sub_id_str = num2str(sub_id_raw);
    elseif isstring(sub_id_raw)
        sub_id_str = char(sub_id_raw);
    elseif ischar(sub_id_raw)
        sub_id_str = sub_id_raw;
    elseif iscategorical(sub_id_raw)
        sub_id_str = char(sub_id_raw);
    else
        sub_id_str = sprintf('unknownID_%04d', i);
    end

    fprintf('\nProcessing lesion %d / %d: %s\n', ...
        i, height(data), sub_id_str);

    % Read the lesion voxel indices from column 3 onward.
    lesion_cells = table2cell(data(i, 3:end));
    lesion_indices = [];

    for k = 1:numel(lesion_cells)
        value = lesion_cells{k};

        if isnumeric(value) && isscalar(value) && ~isnan(value)
            lesion_indices(end + 1, 1) = value; %#ok<SAGROW>
        elseif ischar(value) || (isstring(value) && isscalar(value))
            numeric_value = str2double(value);

            if ~isnan(numeric_value)
                lesion_indices(end + 1, 1) = numeric_value; %#ok<SAGROW>
            end
        end
    end

    if isempty(lesion_indices)
        warning('No voxel indices were found for lesion: %s', sub_id_str);
        continue;
    end

    % Remove duplicate voxel indices while preserving their original order.
    lesion_indices = unique(lesion_indices, 'stable');

    % Map the whole-volume linear indices to rows of the mean FC matrix.
    [is_valid, indices_in_valid] = ismember(lesion_indices, valid_indices);
    indices_in_valid = indices_in_valid(is_valid);

    if isempty(indices_in_valid)
        warning('No valid voxels were found for lesion: %s', sub_id_str);
        continue;
    end

    n_excluded = numel(lesion_indices) - numel(indices_in_valid);
    if n_excluded > 0
        warning('%d voxel(s) from lesion %s were outside valid_indices.', ...
            n_excluded, sub_id_str);
    end

    % Extract and average the FC maps of all valid lesion voxels.
    lesion_fc_matrix = rho_z_mean(indices_in_valid, :);
    lesion_r_map = mean(lesion_fc_matrix, 1, 'omitnan');

    % Replace non-finite values before saving.
    lesion_r_map(~isfinite(lesion_r_map)) = 0;

    output_file = fullfile(output_dir, ...
        ['Amnesia_rmap_meanFC_', sub_id_str, '.mat']);

    save(output_file, 'lesion_r_map');

    fprintf('Saved %d valid lesion voxels to:\n%s\n', ...
        numel(indices_in_valid), output_file);
    toc;
end

fprintf('\nAll lesions have been processed.\n');
