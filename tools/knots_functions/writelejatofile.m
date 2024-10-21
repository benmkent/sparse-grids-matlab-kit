% Script to write Leja points and quadrature weights to a file xw.txt
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

N_max = 100;
a=-1;
b=1;

% Open the file for writing
fileID = fopen('xw.txt', 'w');

% Write the switch-case structure
fprintf(fileID, 'switch n\n');

% Loop over n from 1 to N_max
for n = 1:N_max
    % Get the precomputed points and weights
    [leja_points, weights] = knots_leja_generator(n, a, b);

    % Write the case for n
    fprintf(fileID, '    case %d\n', n);
    fprintf(fileID, '        x = [');
    fprintf(fileID, '%.16f ;', leja_points);
    fprintf(fileID, '];\n');
    fprintf(fileID, '        w = [');
    fprintf(fileID, '%.16f ', weights);
    fprintf(fileID, '];\n');
end

% Write the end of the switch statement
fprintf(fileID, '    otherwise\n');
fprintf(fileID, '        error(''n must be between 1 and %d'');\n', N_max);
fprintf(fileID, 'end\n');

% Close the file
fclose(fileID);

disp('File "xw.txt" written successfully.');
