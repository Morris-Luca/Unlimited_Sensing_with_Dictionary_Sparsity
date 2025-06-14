function create_plot(directory_path, layout)
% CREATE_PLOT - Generates a combined plot from .mat files in the specified directory.
%
% INPUT:
%   directory_path - Path to the directory containing .mat files.
%   layout         - Layout option for the plot:
%                    'portrait'  - 5 rows, 3 columns (default)
%                    'landscape' - 3 rows, 6 columns
%
% OUTPUT:
%   Saves the combined plot as 'combined_plot.png' in a 'plots' folder within the directory.

    % Validate layout input and set default to 'portrait' if not provided
    if nargin < 2
        layout = 'portrait'; % Default layout
    end
    % Ensure the layout is either 'portrait' or 'landscape'
    if ~ismember(layout, {'portrait', 'landscape'})
        error('Invalid layout option. Choose either "portrait" or "landscape".');
    end

    % Get all .mat files in the specified directory
    mat_files = dir(fullfile(directory_path, '*.mat'));
    % Ensure exactly 6 .mat files are present in the directory
    if length(mat_files) ~= 6
        error('The directory must contain exactly 6 .mat files. Found %d.', length(mat_files));
    end

    % Create a folder for saving plots if it doesn't already exist
    plots_folder = fullfile(directory_path, 'plots');
    if ~exist(plots_folder, 'dir')
        mkdir(plots_folder);
    end

    % Load data from each .mat file into a cell array
    data_structs = cell(1, 6); % Preallocate cell array for data
    file_names = strings(1, 6); % Preallocate array for file names
    for i = 1:6
        mat_file_path = fullfile(mat_files(i).folder, mat_files(i).name); % Full path to the .mat file
        data_structs{i} = load(mat_file_path); % Load the .mat file
        file_names(i) = mat_files(i).name; % Store the file name
    end

    % Set up the figure and layout based on the chosen format
    figure;
    if strcmp(layout, 'portrait')
        % Portrait layout: 5 rows, 3 columns
        set(gcf, 'Position', [100, 100, 750, 2000]); % Tall and narrow figure
        t = tiledlayout(5, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    else
        % Landscape layout: 3 rows, 5 columns
        set(gcf, 'Position', [100, 100, 2000, 750]); % Wide and short figure
        t = tiledlayout(3, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
    end

    % Define plot titles for each matrix type
    title = {'Total sparsity', 'Rel. err OMP', 'Rel. err BP', 'Error >0.05 OMP', 'Error >0.05 BP', ''};

    % Loop through columns (1 to 3) to generate plots
    for col = 1:3
        title_idx = 6; % Default title index
        row = 0; % Initialize row index for each column
        for matrix_idx = 1:3 % Loop through matrix types (sparsity, error, success)

            % Define color axis ranges for each matrix type
            caxis_range = cell(1, 3); % Preallocate cell array for color axis ranges
            if col == 1
                caxis_range{1} = [0 400]; % Special range for sparsity matrix in column 1
            else
                caxis_range{1} = [0 100]; % Default range for sparsity matrix
            end
            caxis_range{2} = [0 0.4]; % Range for error matrix
            caxis_range{3} = [0 1]; % Range for success matrix

            for i = 1:2 % Loop through two iterations (e.g., OMP and BP)
                % Skip the second iteration for sparsity matrix
                if matrix_idx == 1 && i == 2
                    continue
                end

                % Load data for the current column and iteration
                data = data_structs{col + 3 * (i - 1)}; % Select the appropriate data structure

                % Extract axes and matrices from the loaded data
                sparsity_axis = data.sparsity_axis; % Sparsity levels
                folding_axis = data.folding_axis; % Folding heights
                matrix_structs = cell(1, 3); % Preallocate cell array for matrices
                matrix_structs{1} = data.sparsity_matrix; % Sparsity matrix
                matrix_structs{2} = data.error_matrix; % Error matrix
                matrix_structs{3} = 1 - data.success_matrix; % Success matrix (inverted)

                % Determine the title and plot position based on layout
                if strcmp(layout, 'portrait')
                    if col == 2
                        title_idx = row + 1; % Update title index for column 2
                    end
                    % Plot the matrix in the appropriate tile
                    plot_matrix(t, 3 * row + col, sparsity_axis, folding_axis, matrix_structs{matrix_idx}, title{title_idx}, caxis_range{matrix_idx}, layout);
                elseif strcmp(layout, 'landscape')
                    if col == 1
                        title_idx = row + 1; % Update title index for column 1
                    end
                    % Plot the matrix in the appropriate tile
                    plot_matrix(t, 5 * (col - 1) + row + 1, sparsity_axis, folding_axis, matrix_structs{matrix_idx}, title{title_idx}, caxis_range{matrix_idx}, layout);
                end
                row = row + 1; % Increment row index
            end
        end
    end

    % Save the combined plot as a PNG file
    saveas(gcf, fullfile(plots_folder, 'combined_plot.png'));
    close(gcf); % Close the figure
end

function plot_matrix(t, tile_index, x_axis, y_axis, matrix_data, plot_title, caxis_range, layout)
% PLOT_MATRIX - Helper function to plot a single matrix in a tiled layout.
%
% INPUT:
%   t            - Tiled layout object.
%   tile_index   - Index of the tile in the layout.
%   x_axis       - X-axis labels (e.g., sparsity levels).
%   y_axis       - Y-axis labels (e.g., folding heights).
%   matrix_data  - Matrix to be visualized.
%   plot_title   - Title of the plot.
%   caxis_range  - Color axis range for the plot.
%   layout       - Layout type ('portrait' or 'landscape').

    % Create a new tile in the layout
    ax = nexttile(t, tile_index);
    imagesc(matrix_data); % Display the matrix as an image

    % Set colormap and colorbar
    colormap('parula');
    colorbar;

    % Set color axis range if provided
    if ~isempty(caxis_range)
        caxis(caxis_range);
    end

    % Add title and axis labels
    title(plot_title, 'FontSize', 15);
    set(gca, 'FontSize', 12, 'YDir', 'normal'); % Set font size and axis direction

    % Set x-axis and y-axis ticks and labels
    xticks(linspace(1, length(x_axis), min(5, length(x_axis))));
    yticks(linspace(1, length(y_axis), min(5, length(y_axis))));
    xticklabels(x_axis(round(linspace(1, length(x_axis), min(5, length(x_axis))))));
    yticklabels(y_axis(round(linspace(1, length(y_axis), min(5, length(y_axis))))));

    % Add axis labels for specific tiles
    if strcmp(layout, 'portrait')
        if mod(tile_index, 3) == 1
            ylabel('Folding height', 'FontSize', 12);
        end
        if tile_index > 12
            xlabel('Sparsity level', 'FontSize', 12);
        end
    elseif strcmp(layout, 'landscape')
        if mod(tile_index, 5) == 1
            ylabel('Folding height', 'FontSize', 12);
        end
        if tile_index > 10
            xlabel('Sparsity level', 'FontSize', 12);
        end
    end

    % Force square axis
    axis square;
end