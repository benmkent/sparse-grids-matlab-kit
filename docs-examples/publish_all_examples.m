filePath = 'sparse_grids_tutorial.m';
options.format = 'html';                 % Output format as HTML
options.outputDir = ('../docs'); % Directory for published files
options.figureSnapMethod = 'print';      % Use the print method for figures
options.maxOutputLines = 10;             % Limit output to 10 lines
publish(filePath, options);

filePath = 'test_compute_normal_leja_and_convergence_test.m';
options.format = 'html';                 % Output format as HTML
options.outputDir = ('../docs'); % Directory for published files
options.figureSnapMethod = 'print';      % Use the print method for figures
options.maxOutputLines = 10;             % Limit output to 10 lines
publish(filePath, options);

filePath = 'test_convert_to_modal.m';
options.format = 'html';                 % Output format as HTML
options.outputDir = ('../docs'); % Directory for published files
options.figureSnapMethod = 'print';      % Use the print method for figures
options.maxOutputLines = 10;             % Limit output to 10 lines
publish(filePath, options);

filePath = 'test_evaluate_on_sparse_grids';
options.format = 'html';                 % Output format as HTML
options.outputDir = ('../docs'); % Directory for published files
options.figureSnapMethod = 'print';      % Use the print method for figures
options.maxOutputLines = 10;             % Limit output to 10 lines
publish(filePath, options);

filePath = 'test_sparse_interpolation';
options.format = 'html';                 % Output format as HTML
options.outputDir = ('../docs'); % Directory for published files
options.figureSnapMethod = 'print';      % Use the print method for figures
options.maxOutputLines = 10;             % Limit output to 10 lines
publish(filePath, options);

filePath = 'test_sparse_quadrature';
options.format = 'html';                 % Output format as HTML
options.outputDir = ('../docs'); % Directory for published files
options.figureSnapMethod = 'print';      % Use the print method for figures
options.maxOutputLines = 10;             % Limit output to 10 lines
publish(filePath, options);

filePath = 'test_spectral';
options.format = 'html';                 % Output format as HTML
options.outputDir = ('../docs'); % Directory for published files
options.figureSnapMethod = 'print';      % Use the print method for figures
options.maxOutputLines = 10;             % Limit output to 10 lines
publish(filePath, options);

filePath = 'tutorial_adaptive';
options.format = 'html';                 % Output format as HTML
options.outputDir = ('../docs'); % Directory for published files
options.figureSnapMethod = 'print';      % Use the print method for figures
options.maxOutputLines = 10;             % Limit output to 10 lines
publish(filePath, options);

filePath = 'tutorial_PlateauSC_adaptive';
options.format = 'html';                 % Output format as HTML
options.outputDir = ('../docs'); % Directory for published files
options.figureSnapMethod = 'print';      % Use the print method for figures
options.maxOutputLines = 10;             % Limit output to 10 lines
publish(filePath, options);

% Specify the directory containing the .html files
folder = '../docs/';
    
% Get a list of all .html files in the folder
files = dir(folder, '*.html');

% Loop through each .html file and copy it with a .md extension
for i = 1:length(files)
    % Get the full file path for the .html file
    htmlFile = fullfile(folder, files(i).name);
    
    % Construct the new file path with .md extension
    [~, name, ~] = fileparts(files(i).name);
    mdFile = fullfile(folder, [name, '.md']);
    
    % Copy the content of the .html file to the .md file
    copyfile(htmlFile, mdFile);
    
    fprintf('Copied: %s to %s\n', htmlFile, mdFile);
end