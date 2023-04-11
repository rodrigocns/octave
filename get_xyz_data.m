#returns cell array with element symbols and xyz coordinates matrix(nx3) of atoms in .xyz file
function [atom_count,elem,atom_coords] = get_xyz_data (filename) 
  fid = fopen(filename, 'r');

  # Read the number of atoms from the first line of the file
  line = fgetl(fid);
  atom_count = str2num(line);

  # Pre-allocate arrays for the element symbols and coordinates
  elem = cell(atom_count, 1);
  atom_coords = zeros(atom_count, 3);

  # Loop over each line of the file after the first line
  line = fgetl(fid); #
  for i = 1:atom_count
    line = fgetl(fid);
    line_data = strsplit(line);
    elem{i} = line_data{1};
    atom_coords(i, :) = str2double(line_data(2:4));
  end

  fclose(fid);

endfunction