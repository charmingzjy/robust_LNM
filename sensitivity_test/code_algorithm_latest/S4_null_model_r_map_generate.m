clear; clc;

%% ===== Parameters =====
N_NULL = 10000;
START_NULL = 1;  % Change this value to resume an interrupted run

% Select the S3 null-model source:
%   'annealing'  - frequency-matched null models
%   'reposition' - randomly repositioned null models
NULL_MODEL_METHOD = 'annealing';

% When true, existing lesion r-map files are retained and skipped.
SKIP_EXISTING_FILES = true;

%% ===== Project paths =====
script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(script_dir);

template_dir = fullfile(project_dir, 'template');
project_output_dir = fullfile(project_dir, 'output');

valid_indices_file = fullfile( ...
    template_dir, 'valid_indices_standard_wh_MNI8mm_32+_32_32.mat');
mean_fc_file = fullfile(template_dir, 'rho_z_mean.mat');

switch lower(NULL_MODEL_METHOD)
    case 'annealing'
        null_input_dir = fullfile( ...
            project_output_dir, ...
            sprintf( ...
            'R2_null_model_random_indices_iter%d_sortedfreq', N_NULL));
        output_root = fullfile( ...
            project_output_dir, ...
            sprintf('R3_null%d_generate_r_map_sortedfreq', N_NULL));

    case 'reposition'
        null_input_dir = fullfile( ...
            project_output_dir, ...
            sprintf( ...
            'R2_null_model_random_indices_iter%d_reposition', N_NULL));
        output_root = fullfile( ...
            project_output_dir, ...
            sprintf('R3_null%d_generate_r_map_reposition', N_NULL));

    otherwise
        error(['NULL_MODEL_METHOD must be either "annealing" ' ...
            'or "reposition".']);
end

%% ===== Validate parameters and input paths =====
if START_NULL < 1 || START_NULL ~= floor(START_NULL)
    error('START_NULL must be a positive integer.');
end

if N_NULL < START_NULL || N_NULL ~= floor(N_NULL)
    error('N_NULL must be an integer greater than or equal to START_NULL.');
end

if ~exist(null_input_dir, 'dir')
    error('S3 null-model directory not found: %s', null_input_dir);
end

if ~isfile(valid_indices_file)
    error('Valid-indices file not found: %s', valid_indices_file);
end

if ~isfile(mean_fc_file)
    error('Mean FC file not found: %s', mean_fc_file);
end

if ~exist(output_root, 'dir')
    mkdir(output_root);
end

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

if any(~isfinite(valid_indices) | valid_indices < 1 | ...
        valid_indices ~= floor(valid_indices))
    error('Variable "valid_indices" contains invalid voxel indices.');
end

if numel(unique(valid_indices)) ~= numel(valid_indices)
    error('Variable "valid_indices" contains duplicate voxel indices.');
end

fprintf('Loaded %d valid voxel indices.\n', numel(valid_indices));

%% ===== Load the GSP1000 mean functional-connectivity matrix =====
mean_fc_data = load(mean_fc_file);

if isfield(mean_fc_data, 'rho_z')
    rho_z_mean = mean_fc_data.rho_z;
elseif isfield(mean_fc_data, 'rho_z_mean')
    rho_z_mean = mean_fc_data.rho_z_mean;
else
    error('Neither "rho_z" nor "rho_z_mean" was found in: %s', ...
        mean_fc_file);
end

if ~isnumeric(rho_z_mean) || isempty(rho_z_mean) || ...
        ~ismatrix(rho_z_mean)
    error('The mean FC matrix must be a nonempty numeric matrix.');
end

if size(rho_z_mean, 1) ~= numel(valid_indices)
    error(['The number of mean FC rows (%d) does not match the ' ...
        'number of valid indices (%d).'], ...
        size(rho_z_mean, 1), numel(valid_indices));
end

fprintf('Mean FC matrix loaded. Size = %d x %d.\n', ...
    size(rho_z_mean, 1), size(rho_z_mean, 2));

%% ===== Generate r-maps for all null lesions =====
for null_number = START_NULL:N_NULL
    fprintf('\n=== Processing null %05d / %05d ===\n', ...
        null_number, N_NULL);

    null_mat_file = fullfile( ...
        null_input_dir, sprintf('null_%05d.mat', null_number));

    if ~isfile(null_mat_file)
        warning('Null-model file not found; skipped: %s', null_mat_file);
        continue;
    end

    null_file_data = load(null_mat_file, 'null_data');

    if ~isfield(null_file_data, 'null_data') || ...
            ~isstruct(null_file_data.null_data)
        warning('Variable "null_data" was not found or invalid: %s', ...
            null_mat_file);
        continue;
    end

    null_data = null_file_data.null_data;

    output_folder = fullfile( ...
        output_root, sprintf('null_%05d', null_number));

    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    num_saved = 0;
    num_skipped = 0;

    for i = 1:numel(null_data)
        if ~isfield(null_data, 'Subid') || ...
                ~isfield(null_data, 'lin_idx')
            error('Required fields are missing from null_data in: %s', ...
                null_mat_file);
        end

        subject_id = value_to_text(null_data(i).Subid, i);
        safe_subject_id = make_safe_filename(subject_id);

        output_name = sprintf( ...
            'Amnesia_rmap_meanFC_null%05d_%s.mat', ...
            null_number, safe_subject_id);
        output_file = fullfile(output_folder, output_name);

        if SKIP_EXISTING_FILES && isfile(output_file)
            num_skipped = num_skipped + 1;
            continue;
        end

        lesion_indices = double(null_data(i).lin_idx(:));

        if isempty(lesion_indices)
            warning('Null %05d, lesion %s has no voxel indices.', ...
                null_number, subject_id);
            continue;
        end

        if any(~isfinite(lesion_indices) | lesion_indices < 1 | ...
                lesion_indices ~= floor(lesion_indices))
            warning('Null %05d, lesion %s has invalid voxel indices.', ...
                null_number, subject_id);
            continue;
        end

        % Map whole-volume linear indices to rows of the mean FC matrix.
        [is_valid, indices_in_valid] = ...
            ismember(lesion_indices, valid_indices);
        indices_in_valid = indices_in_valid(is_valid);

        if isempty(indices_in_valid)
            warning('No valid voxels for lesion %s in null %05d.', ...
                subject_id, null_number);
            continue;
        end

        num_excluded = numel(lesion_indices) - numel(indices_in_valid);
        if num_excluded > 0
            warning(['Null %05d, lesion %s: %d voxel(s) were outside ' ...
                'valid_indices.'], ...
                null_number, subject_id, num_excluded);
        end

        lesion_fc_matrix = rho_z_mean(indices_in_valid, :);
        lesion_r_map = mean(lesion_fc_matrix, 1, 'omitnan');
        lesion_r_map(~isfinite(lesion_r_map)) = 0;

        save(output_file, 'lesion_r_map');
        num_saved = num_saved + 1;
    end

    fprintf('Null %05d complete: %d saved, %d existing skipped.\n', ...
        null_number, num_saved, num_skipped);
end

fprintf('\nAll available null-model r-maps were processed successfully.\n');

%% ===== Local functions =====
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

function safe_text = make_safe_filename(text)
    safe_text = regexprep(strtrim(text), '[<>:"/\\|?*]', '_');
    safe_text = regexprep(safe_text, '\s+', '_');

    if isempty(safe_text)
        safe_text = 'unknownID';
    end
end
