function speciesStruct = packSpeciesStruct(speciesNames, values)
    % FUNCTION NAME:
    %   unpackIntMatrix
    %
    % DESCRIPTION:
    %   Create variable names for the entries in the interference matrix in
    %   the parameters container map
    %
    % INPUT:
    %   speciesNames - (string array) All of the species names (in the
    %       order of the values stored in x)
    %   values (double array) The values associated with the different
    %       species (x in the ODE function)
    %
    % OUTPUT:
    %   speciesStruct - (struct) The species names and their value
    %
    % ASSUMPTIONS AND LIMITATIONS:
    %
    % REVISION HISTORY:
    %   2023-09-05 - Helen Scott
    %       * Initial implementation
    
    speciesStruct = struct();
    for idx = 1:length(speciesNames)
        speciesStruct.(speciesNames(idx)) = values(idx);
    end