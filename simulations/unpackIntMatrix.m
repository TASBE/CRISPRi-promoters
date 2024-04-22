function parameters = unpackIntMatrix(parameters, i_matrix)
    % FUNCTION NAME:
    %   unpackIntMatrix
    %
    % DESCRIPTION:
    %   Create variable names for the entries in the interference matrix in
    %   the parameters container map
    %
    % INPUT:
    %   parameters - (map) Map object listing the parameter names and
    %       values as determined from parameter fitting and literature
    %   i_matrix - (double) 2D interference matrix 
    %
    % OUTPUT:
    %   parameters - (map) Map object listing the parameter names and
    %       values as determined from parameter fitting and literature,
    %       including all i_matrix values
    %
    % ASSUMPTIONS AND LIMITATIONS:
    %
    % REVISION HISTORY:
    %   2023-01-31 - Helen Scott
    %       * Initial implementation

    for row = 1:size(i_matrix, 1)
        for col = 1:size(i_matrix, 2)
             parameters("int_matrix_" + int2str(row) + "_" + int2str(col)) = i_matrix(row, col);
        end
    end