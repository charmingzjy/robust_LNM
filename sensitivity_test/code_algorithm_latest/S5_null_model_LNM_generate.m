clear; clc;

%% ===== Parameters =====
N_NULL = 10000;
START_NULL = 1;  % Change this value to resume an interrupted run

% Select the S4 r-map source:
%   'annealing'  - frequency-matched null models
%   'reposition' - randomly repositioned null models
NULL_MODEL_METHOD = 'annealing';

% These thresholds must match S2_true_lesion_LNM_generate.m.
LOWER_PERCENTILE = 30;
UPPER_PERCENTILE = 70;

% When true, existing null LNM result files are retained and skipped.
SKIP_EXISTING_FILES = true;

%% ===== Project paths =====

script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(script_dir);
project_output_dir = fullfile(project_dir, 'output');

switch lower(NULL_MODEL_METHOD)
    case 'annealing'
        null_rmap_root = fullfile( ...
            project_output_dir, ...
            sprintf('R3_null%d_generate_r_map_sortedfreq', N_NULL));
        output_dir = fullfile( ...
            project_output_dir, ...
            sprintf(['R4_LNM_rmaps_null_%d_results_' ...
            'percent30_sortedfreq'], N_NULL));

    case 'reposition'
        null_rmap_root = fullfile( ...
            project_output_dir, ...
            sprintf('R3_null%d_generate_r_map_reposition', N_NULL));
        output_dir = fullfile( ...
            project_output_dir, ...
            sprintf(['R4_LNM_rmaps_null_%d_results_' ...
            'percent30_reposition'], N_NULL));

    otherwise
        error(['NULL_MODEL_METHOD must be either "annealing" ' ...
            'or "reposition".']);
end

%% ===== Validate parameters and paths =====
if START_NULL < 1 || START_NULL ~= floor(START_NULL)
    error('START_NULL must be a positive integer.');
end

if N_NULL < START_NULL || N_NULL ~= floor(N_NULL)
    error('N_NULL must be an integer greater than or equal to START_NULL.');
end

if LOWER_PERCENTILE < 0 || UPPER_PERCENTILE > 100 || ...
        LOWER_PERCENTILE >= UPPER_PERCENTILE
    error('The percentile thresholds are invalid.');
end

if ~exist(null_rmap_root, 'dir')
    error('S4 null r-map directory not found: %s', null_rmap_root);
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% ===== Generate one signed overlap map for each null dataset =====
for null_number = START_NULL:N_NULL
    output_file = fullfile( ...
        output_dir, ...
        sprintf('Amnesia_LNM_rmap_null%05d.mat', null_number));

    if SKIP_EXISTING_FILES && isfile(output_file)
        fprintf('Null %05d already exists; skipped.\n', null_number);
        continue;
    end

    null_folder = fullfile( ...
        null_rmap_root, sprintf('null_%05d', null_number));

    if ~exist(null_folder, 'dir')
        warning('Null r-map folder not found; skipped: %s', null_folder);
        continue;
    end

    file_pattern = sprintf( ...
        'Amnesia_rmap_meanFC_null%05d_*.mat', null_number);
    files = dir(fullfile(null_folder, file_pattern));

    if isempty(files)
        % This fallback also accepts older S4 output naming.
        files = dir(fullfile( ...
            null_folder, 'Amnesia_rmap_meanFC*.mat'));
    end

    num_files = numel(files);

    if num_files == 0
        warning('No lesion r-map files found; skipped: %s', null_folder);
        continue;
    end

    % Use a consistent processing order.
    [~, sort_order] = sort({files.name});
    files = files(sort_order);

    fprintf('\n=== Processing null %05d / %05d (%d lesions) ===\n', ...
        null_number, N_NULL, num_files);

    %% ===== Load all lesion r-maps for this null dataset =====
    first_r_map = load_r_map(fullfile(null_folder, files(1).name));
    num_voxels = numel(first_r_map);

    r_maps = zeros(num_voxels, num_files);
    r_maps(:, 1) = first_r_map;

    for j = 2:num_files
        file_path = fullfile(null_folder, files(j).name);
        r_map = load_r_map(file_path);

        if numel(r_map) ~= num_voxels
            error(['R-map length mismatch in "%s": expected %d values, ' ...
                'but found %d.'], ...
                files(j).name, num_voxels, numel(r_map));
        end

        r_maps(:, j) = r_map;
    end

    %% ===== Select the top and bottom percent for each lesion =====
    masks = zeros(num_voxels, num_files, 'int8');

    for j = 1:num_files
        current_r_map = r_maps(:, j);

        upper_threshold = prctile( ...
            current_r_map, UPPER_PERCENTILE);
        lower_threshold = prctile( ...
            current_r_map, LOWER_PERCENTILE);

        masks(current_r_map >= upper_threshold, j) = 1;
        masks(current_r_map <= lower_threshold, j) = -1;
    end

    %% ===== Determine the signed majority at each voxel =====
    positive_count = sum(masks == 1, 2);
    negative_count = sum(masks == -1, 2);
    majority_difference = positive_count - negative_count;

    majority_count = zeros(num_voxels, 1);

    positive_majority = majority_difference > 0;
    negative_majority = majority_difference < 0;

    majority_count(positive_majority) = ...
        positive_count(positive_majority);
    majority_count(negative_majority) = ...
        -negative_count(negative_majority);

    %% ===== Calculate and save the signed overlap ratio =====
    % Positive values indicate a positive-connectivity majority.
    % Negative values indicate a negative-connectivity majority.
    % Zero indicates a tie or no majority.
    overlap = majority_count / num_files;

    save(output_file, 'overlap');

    fprintf('Null %05d saved to:\n%s\n', null_number, output_file);
end

fprintf('\nAll available null datasets were processed successfully.\n');

%% ===== Local function =====
function r_map = load_r_map(file_path)
    file_data = load(file_path);

    if isfield(file_data, 'lesion_r_map')
        r_map = file_data.lesion_r_map;
    elseif isfield(file_data, 'r_values')
        r_map = file_data.r_values;
    else
        error('Neither "lesion_r_map" nor "r_values" was found in: %s', ...
            file_path);
    end

    if ~isnumeric(r_map) || isempty(r_map) || ~isvector(r_map)
        error('The r-map must be a nonempty numeric vector: %s', ...
            file_path);
    end

    r_map = double(r_map(:));

    if any(~isfinite(r_map))
        warning('Non-finite values in "%s" were replaced with zero.', ...
            file_path);
        r_map(~isfinite(r_map)) = 0;
    end
end
