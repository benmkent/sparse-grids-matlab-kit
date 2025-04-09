
%% Publish
addpath(genpath('../'))

close all
filePath = 'sparse_grids_tutorial.m';
publish(filePath, 'format','html','stylesheet','stylesheet.xsl','maxOutputLines',20,'figureSnapMethod','print','outputDir','../docs','catchError',true);

close all
filePath = 'test_spectral';
publish(filePath, 'format','html','stylesheet','stylesheet.xsl','maxOutputLines',20,'figureSnapMethod','print','outputDir','../docs','catchError',true);

close all
filePath = 'tutorial_adaptive';
publish(filePath, 'format','html','stylesheet','stylesheet.xsl','maxOutputLines',20,'figureSnapMethod','print','outputDir','../docs','catchError',true);

close all
filePath = 'tutorial_PlateauSC_adaptive';
publish(filePath, 'format','html','stylesheet','stylesheet.xsl','maxOutputLines',20,'figureSnapMethod','print','outputDir','../docs','catchError',true);

%% Copy to .md

% Specify the directory containing the .html files
folder = '../docs/';
    
% Get a list of all .html files in the folder
files = dir(fullfile(folder, '*.html'));

% Loop through each .html file and copy it with a .md extension
for i = 1:length(files)
    % Get the full file path for the .html file
    htmlFile = fullfile(folder, files(i).name);
    
    % Construct the new file path with .md extension
    [~, name, ~] = fileparts(files(i).name);
    txtFile = fullfile(folder, [name, '.txt']);
    mdFile = fullfile(folder, [name, '.md']);
    
    % Copy the content of the .html file to the .md file
    copyfile(htmlFile, txtFile);

    fid = fopen(txtFile, 'rt');  % Open the file for reading (text mode)
    if fid == -1
        error('File not found or cannot be opened.');
    end

    fileContent = fread(fid, '*char')';  % Read the file contents as a string
    fclose(fid);  % Close the file
    % Perform find and replace
    oldString = 'ENDHTML ';  % The text you want to find
    newString = ['</html>',newline];  % The text you want to replace it with
    fileContent = strrep(fileContent, oldString, newString);  % Replace text
    oldString = 'STARTHTML';  % The text you want to find
    newString = [newline,'<html>'];  % The text you want to replace it with
    fileContent = strrep(fileContent, oldString, newString);  % Replace text

    % Write the modified content back to the file
    fid = fopen(txtFile, 'wt');  % Open the file for writing (text mode)
    if fid == -1
        error('File cannot be opened for writing.');
    end

    fprintf(fid, '%s', fileContent);  % Use %s to write the string with the newline
    fclose(fid);  % Close the file
    copyfile(txtFile, mdFile);
    delete(txtFile);
    delete(htmlFile);

    fprintf('Copied: %s to %s\n', htmlFile, mdFile);
end