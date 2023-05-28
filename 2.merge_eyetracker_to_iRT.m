pkg load io
%{
This script atempts to read the eyeTracker and iRT files to pre-process them
and then merge both data in time. We used the UNIX EPOCH time to do this, as it
is the innate time used in most programming languages, especially javascript and
HTML. If the data from your eyetracker device/software does not contain this
kind of information, then it must be added using the epoch data produced from
iRT, or changing from another epoch system into the unix epoch.

Maybe usefull links:
https://wiki.octave.org/IO_package
%}

% clear all but the slowest variables to obtain
clear -exclusive *_data;

% SETTINGS
% leave 'true' to let the script calculate (or recalculate/read/write again) the described values

% read eyeTracking .xlsx input file
cfg_eyeT_input = false; %slow
cfg_eyeT_input_filename = "raw_eyeT_r.xlsx"; %name of the input file (.xlsx, numbers only, no commas for decimals)

% (only if necessary) append a column of epoch unix data to the eyeT data array
cfg_add_epoch_from_milliseconds = true;
cfg_epoch_ms_col = 3; %time column in milliseconds used to create epoch column
cfg_epoch_anchor = 1682707674981-5656; %estimated unix epoch at time 0 inside eyetracking data.

% fix missing pupil data by linear interpolation (best to always leave on with a new data arrays)
cfg_pupil_fix = true;
cfg_pupil_cols = [9,10]; % eyeT_data columns that need the interpolation. Pulil diameter, left/right, etc.
cfg_pupil_missing_data = [0,-1]; % possible results of missing data to be identified and corrected

% read iRT .xlsx input file
cfg_iRT_input = false; %slow
cfg_iRT_input_filename = "iRT_data.xlsx"; %name of the iRT sheets file unpackaged by unpacking_sheets.m

% merge eyeTracker data to the chosen iRT task (be sure that the files are from the same session)
cfg_data_merge = true;
cfg_iRT_sessionID = 1682707472090; %session ID of the desired task
cfg_iRT_taskID = "bolaBastao_c"; %task ID of the desired task
cfg_iRT_cols = [3:8]; %range of desired data columns from raw_iRT_data. 5:8 is quaternion data, 3 is unix epoch
cfg_eyeT_cols = [1,2,4,7:10]; %range of desired data columns from eyeT_data. epoch data was appended as the first column
cfg_plot_resolugram = false; % DRAW resolugram plot (distance between reference and interactive models, in degrees) (figure#1)

% WRITE output file from task_data (BUGged)
cfg_write_output = false;
cfg_output_filename = strcat("mergeOutput_", cfg_iRT_taskID, num2str(cfg_iRT_sessionID), ".xlsx"); % name of the output file (.xlsx) from the merging of eyeT_data and iRT_data

% read .xyz file with atom data from the used model based on values in session_data
cfg_xyz_input = true;
cfg_xyz_col = 11; % index of the column to look for the modelName value in session_data
cfg_xyz_plot = false; % DRAW scatter3 of the array of atoms colored acording to atom_elem (figure#2)

% calculate temporal array of rotated atoms
cfg_atom_matrix = true;
cfg_atom_matrix_quat_cols = [3:6]; % quaternion index of columns inside task_data array
cfg_atom_matrix_ref_cols = [7:10]; % column indexes for the reference quaternion values inside session_data ("ref_i","ref_j","ref_k","ref_theta")
cfg_atom_matrix_ref_plot = false; %2d scatter of the reference model (figure#3)
cfg_atom_matrix_plot = false; %2d scatter of the interactive model at a set time (figure#4)
cfg_atom_matrix_plot_t = 1; %frame used in rotated cfg_atom_matrix_plot (figure#4)

% calculate gaze related things
cfg_gaze_dist_matrix = true;
cfg_gaze_cols = [10,11]; % column indexes of gaze x and y coordinates on screen. Should be in pixels, counting from top-left corner
cfg_gaze_scrSize_cols = [12,13]; % column indexes of screenSize values from session_data (width and height respectively)
cfg_gaze_cvsRef_cols = [14,15,16,17]; % column indexes of reference model canvas positionsfrom session_data (in order: top, right, bottom, left)
cfg_gaze_cvsInt_cols = [18,19,20,21]; % column indexes of interactive model canvas positions from session_data (in order: top, right, bottom, left)
cfg_gaze_pxAngs_rate_col = [6]; %column index of pixels (screen distance) per angstrom (atomic distance unit in jmol) in session_data.

%cfg_gaze_
cfg_gaussian_wdt = 0.9; % width of gaussian
cfg_gaze_status_codeInt = 2; %condition in gaze_status, meaning that gaze was within Interactive model canvas
cfg_gaze_status_codeRef = 1; %condition in gaze_status, meaning that gaze was within Reference model canvas
cfg_gaze_status_array = true;

cfg_gaze_heatmap = true;
cfg_gaze_heatmap_window = true;
cfg_heatmap_mw_frame_length = 20; %moving window length in frames for heatmap computation
cfg_gaussian_wdt = 150; %gaussian width in screen pixels, used in heatmap calculation
cfg_heatmap_mw_filename = "translucent temporal heatmap.xlsx";
cfg_heatmap_filename = "translucent heatmap.xlsx";

%{
   #=========================================#
   # DON'T MODIFY ANYTHING BELLOW THIS LINE! #
   #     (maybe add_epoch() function if)     #
   #      (you know what you are doing)      #
   #=========================================#
%}

%FUNCTIONS
% recognize input file format and reads eyeTracker data (eyeT_data)
function [raw_eyeT_data, raw_eyeT_header_data] = eyeT2oct (filename)
  tic();
  printf("Eyetracker data input: opening..");
  %raw_eyeT_data = xlsread (filename);
  xls_eyeT = xlsopen(filename);
  printf(".reading (this can take time, aprox. 1 min for 10k lines,8 cols).");
  cell_data = xls2oct(xls_eyeT);
  printf(".parsing.");
  raw_eyeT_header_data = cell(cell_data(1,:));
  raw_eyeT_data = parsecell(cell_data);
  printf(".closing.");
  [xls_eyeT] = xlsclose(xls_eyeT);
  printf(".done! ");
  toc();
endfunction
% add epoch data to a time column counting milliseconds
function epoch_column = add_epoch (raw_eyeT_data, epoch_column, epoch_anchor)
  printf("Adding epoch column..");
  epoch_column = raw_eyeT_data(:,epoch_column);
  epoch_column = epoch_column + epoch_anchor;
  printf("Done!\n");
endfunction
% interpolate missing eyeT_data
function interpolated_data = interpolate_missing_data (input_data, columns_to_fix, missing_data)
  printf("Interpolating missing values.. ");
  interpolated_data = input_data;
  % loop for each column
  for c_i=1:numel(columns_to_fix)
    col = columns_to_fix(c_i);
    printf(".[col = %i]",col);
    % find first valid value
    first_valid_line = 1;
    while ismember (interpolated_data(first_valid_line,col), missing_data)
      first_valid_line++;
    endwhile
    printf("[1st valid line = %i].",first_valid_line);
    % find missing values beyond 1st valid line/value
    for line=first_valid_line+1:size(input_data,1)
      if ismember (interpolated_data(line,col), missing_data)
        i=1; % amount of sequential missing values
        while ismember (interpolated_data(line+i,col), missing_data)
          i++;
          % IF line+i is outside of data range, it means the last array value is missing
          if line+i > size(interpolated_data,1)
            % make last value equal to previous valid value and continue process
            i--;
            interpolated_data(line+i,col)=interpolated_data(line-1,col);
            break;
          endif
        endwhile
        %filling values
        first_num = interpolated_data(line-1,col); %line before missing value
        last_num = interpolated_data(line+i,col); %line with value
        clear fix_array;
        fix_array(1:2+i,1) = linspace(first_num, last_num, 2+i); %vertical arr w/ interpolated values
        %printf("[%i]",line);
        interpolated_data(line-1:line+i,col) = fix_array;
      endif
    endfor
  endfor
  printf(".values interpolated!\n");
endfunction
% read data from iRT .xlsx; 'session_arr' has session details, 'raw_arr' has the data bulk in a cell matrix
function [session_arr, raw_arr] = iRT2oct (filename)
  tic();
  printf("Reading iRT data (this can take a while)..");
  xls_iRT = xlsopen(filename);
  session_arr = xls2oct(xls_iRT,"sessions");
  raw_arr = xls2oct(xls_iRT,"data");
  [xls_iRT] = xlsclose(xls_iRT);
  printf(".done!");
  toc();
endfunction
% create data array and its header from slice of iRT data (task data) corresponding to session_ID and task_ID.
function [numeric_data, header_data] = slice_task_data (raw_arr, session_ID, task_ID, desired_iRT_columns)
  % 'desired_iRT_columns' is the range of the desired columns from raw_iRT_data (5:8 or 3,4,5:8)
printf("Slicing out task data..");
  % find first line of the slice
  first_line= -1; %defined outside range
  for n = 1 :size(raw_arr,1)
    if and( strcmp(raw_arr{n,1},num2str(session_ID)), strcmp(raw_arr{n,2}, task_ID) )
      first_line = n;
      printf(".first line of slice is %i.", first_line);
      break;
    endif
  endfor
  % find last line of the slice
  last_line= -1;
  for n = first_line : size(raw_arr,1)
    if not( and( strcmp( raw_arr{n,1}, num2str( session_ID) ), strcmp( raw_arr{n,2}, task_ID) ))
      last_line = n - 1;
      printf(".last line of slice is %i.", last_line);
      break;
    endif
  endfor
  % get slice of data
  numeric_data = cell2mat( raw_arr( first_line:last_line, desired_iRT_columns) );
  header_data = cell(raw_arr(1,desired_iRT_columns));
  printf(".data sliced!\n");
endfunction

% merge eyeT_data into iRT_data based on the nearest time values by nearest neighbours method
function mergedMatrix = merge_data (eyeT_data, task_data, desired_eyeT_columns)
  %Find nearest indexes in eyeT_data for each time point in iRT_data
  printf("Calculating nearest indexes.");
  eyeT_time_col = eyeT_data(:, 1); %eyeT epoch
  iRT_time_col = task_data(:, 1); %iRT epoch
  [~, nearest_indices] = min(abs(eyeT_time_col - iRT_time_col'));
  % Merge desired eyeT_data columns into iRT_data
  printf(".merging to task_data.");
  mergedMatrix = [task_data, eyeT_data(nearest_indices, desired_eyeT_columns)];
  printf(".done!\n");
endfunction

% write output file from merged data
function writeOutput (filename, header, task_data)
  tic(); printf("Writing output file..");
  xls_merged = xlsopen (filename, true);
  merged_cell_matrix = vertcat (header, num2cell (task_data));
  [xls_merged] = oct2xls (merged_cell_matrix , xls_merged, "data");
  [xls_merged] = xlsclose ( xls_merged);
  printf(".done!"); toc();
endfunction
% get rowIndex of session/task pair in session_data.
function rowIndex = get_session_row (session_data, session_ID, task_ID)
  for row = 1 :size(session_data,1)
    if and (strcmp (session_data{row,1}, num2str (session_ID)), strcmp (session_data{row,2}, task_ID) )
      rowIndex = row;
      break;
    endif
  endfor
endfunction
% codify atoms by element|| in:(atom_count,atom_xyz,atom_elem) || out:[size,R,G,B]
function atom_cor = generate_color_vector (atom_count, atom_xyz, atom_elem)
  %used in some graph renderings at the end of this script
  for i = 1:atom_count
    switch ( strvcat(atom_elem(i)) )  %strvcat extracts a string from a cell array
      case "H"
        atom_xyz(i,4) = 1; % atomic number
        atom_cor(i,1:4) = [74, 0.75,0.75,0.75]; %[size, R,G,B]
      case "C"
        atom_xyz(i,4) = 12;
        atom_cor(i,1:4) = [154, 0,0,0];
      case "N"
        atom_xyz(i,4) = 14;
        atom_cor(i,1:4) = [140, 0,0,1];
      case "O"
        atom_xyz(i,4) = 16;
        atom_cor(i,1:4) = [146, 1,0,0];
      case "Co"
        atom_xyz(i,4) = 27;
        atom_cor(i,1:4) = [100, 1,0,0];
      otherwise
        error ("ERRO! Novo elemento detectado. Atualizar codigo");
        atom_xyz(i,4) = -1;
        atom_cor(i,1:4) = [20, 0,1,0];
    endswitch
  endfor
endfunction
% normalize xyz coords so its center of rotation (center of jmol boundingbox ) is 0,0,0
function norm_atom_xyz = normalize_jmol_rot_center (atom_xyz)
  % get max and min x,y,z coords from the atom_xyz array.
  % These are the boundingbox extremities in jmol.
  max_xyz = max(atom_xyz(:,1:3));
  min_xyz = min(atom_xyz(:,1:3));
  % normalize position of entire array
  normalization_center = (max_xyz+min_xyz)/2;
  norm_atom_xyz = atom_xyz(:,1:3) - normalization_center;
endfunction
% read .xyz file with atom data about the model used.
function [atom_count, elem, atom_coords] = get_xyz_data (filename)
  printf( strcat("Opening .xyz file: ", filename, "... ") );
  fid = fopen(filename, 'r');
  % Read the number of atoms from the first line of the file
  printf(" reading...");
  line = fgetl(fid);
  atom_count = str2num(line);

  % Pre-allocate arrays for the element symbols and coordinates
  elem = cell(atom_count, 1);
  atom_coords = zeros(atom_count, 3);

  % Loop over each line of the file after the first line
  line = fgetl(fid);
  for i = 1:atom_count
    line = fgetl(fid);
    line_data = strsplit(line);
    elem{i} = line_data{1};
    atom_coords(i, :) = str2double(line_data(2:4));
  end
  fclose(fid);
  printf(" sucess! \n");
endfunction
% return rotation matrix {R} from a quaternion {i,j,k,real}
function R = rot_matrix (qi,qj,qk,qr, s = 1)
  % real component is last because jmol data uses this order: qi,qj,qk,qr.
  R= [1-2*s*(qj^2 + qk^2), 2*s*(qi*qj - qk*qr), 2*s*(qi*qk + qj*qr);
      2*s*(qi*qj + qk*qr), 1-2*s*(qi^2 + qk^2), 2*s*(qj*qk - qi*qr);
      2*s*(qi*qk - qj*qr), 2*s*(qj*qk + qi*qr), 1-2*s*(qi^2 + qj^2)];
endfunction
% rotate atom_xyz atom matrix in time based on the interactive model quaternions from task_data matrix
function atomInt_xyzRot = rotate_atom_xyz (task_data, atom_matrix_quat_cols, atom_xyz)
  % copy quaternions from interactive model
  Q = task_data(:,atom_matrix_quat_cols);
  % create rotation matrix for each frame (rot_vector)
  for t = 1: size (Q,1)
    % rot_vector(1:3,1:3,1) = [0,1,0;-1,0,0;0,0,1];   %DEBUG: this matrix should rotate atoms in 90degrees
    rot_vector(1:3,1:3,t) = rot_matrix (Q(t,1), Q(t,2), Q(t,3), Q(t,4));
    % apply rotation for each atom.
    for a = 1:size(atom_xyz,1)
      % 3x3 * 3x1 = 3x1 {rotation center at 0,0,0}
      % this is the right order. changing it will give reversed results!
      atomInt_xyzRot(a,1:3,t) = (rot_vector(1:3,1:3,t)*atom_xyz(a,1:3)' )' ;
    endfor
  endfor
endfunction
% rotate atom_xyz atom matrix based on the reference model quaternion from session_data
function atomRef_xyzRot = rotate_ref_atom_xyz (session_data, session_row, cfg_atom_matrix_ref_cols, atom_xyz)
  Q_ref = cell2mat ( session_data(session_row,cfg_atom_matrix_ref_cols) );
  for a = 1:size(atom_xyz,1);
    atomRef_xyzRot(a,1:3) = ( rot_matrix (Q_ref(1), Q_ref(2), Q_ref(3), Q_ref(4)) * atom_xyz(a,1:3)' )' ;
  endfor
endfunction
% calculate, in pixels, the canvas center from the browser window
function canvas_center = get_cvs_center (session_data, row_index, cvs_pos_cols)
  canvas = cell2mat (session_data (row_index,cvs_pos_cols) ); %top, right, bottom, left side positions.
  canvas_center = [ (canvas(2) + canvas(4) )/2 , (canvas(1) + canvas(3) )/2 ];
endfunction
% compute pixel distance between gaze matrix and each atom of selected atom matrix
function [gaze_atomInt_dist, gaze_atomRef_dist] = get_gaze_atom_dist (gaze_px, atom_int_px, atom_ref_px)
%    printf("Calculating gaze-atom distance array.."); tic();
    %initialize array
    atom_count = size(atom_int_px,1);
    row_count = size(gaze_px,1);
    gaze_atom_dist = zeros ( row_count, atom_count );
    % compute matrices of screen position distances between gaze and each atom
    for a=1:atom_count %for each atom
      for t=1 : row_count %for each point in time
        %Formula: gaze_atom_dist = sqrt( (x-x')^2 + (y-y')^2 )
%        gaze_atom_dist(t,a) = sqrt ( (atom_px(a,1,t)-gaze_px(t,1))^2 + (atom_px(a,2,t)-gaze_px(t,2))^2 );
        gaze_atomInt_dist(t,a) = sqrt ( (atom_int_px(a,1,t)-gaze_px(t,1))^2 + (atom_int_px(a,2,t)-gaze_px(t,2))^2 );
        gaze_atomRef_dist(t,a) = sqrt ( (atom_ref_px(a,1)-gaze_px(t,1))^2 + (atom_ref_px(a,2)-gaze_px(t,2))^2 );
      endfor
    endfor
%    printf(".array calculated.");
endfunction
% stores nearest atom distance and index.
function [gaze_nearest_atom_distance,gaze_nearest_atom_index] = fill_nearest_atom (gaze_px, gaze_atomSome_dist)
  frame_count = size(gaze_px, 1);
  %stores nearest atom distance and index.
  for t=1 : frame_count
      [gaze_nearest_atom_distance(t), gaze_nearest_atom_index(t)] =  min (gaze_atomSome_dist(t,:));
  endfor
endfunction

% fill gaze_status (if gaze is inside canvas 1 or 2) from gaze pixel position and canvas extremities position in px [top right bottom left]
function gaze_status = fill_gaze_status (gaze_px, canvas_ref, cfg_gaze_status_codeRef, canvas_int, cfg_gaze_status_codeInt)
  frame_count = size(gaze_px, 1);
  %0,0 is top-left corner
  cvsTop_r = canvas_ref(1);
  cvsRight_r = canvas_ref(2);
  cvsBottom_r = canvas_ref(3);
  cvsLeft_r = canvas_ref(4);
  cvsTop_i = canvas_int(1);
  cvsRight_i = canvas_int(2);
  cvsBottom_i = canvas_int(3);
  cvsLeft_i = canvas_int(4);
  %fill gaze_status if gaze is within canvas bounduaries
  for t=1 : frame_count
    % if gaze is inside reference canvas..
    if ( (cvsLeft_r < gaze_px(t,1)) && (gaze_px(t,1) < cvsRight_r) && (cvsTop_r < gaze_px(t,2)) && (gaze_px(t,2) < cvsBottom_r) )
      gaze_status(t,1) = cfg_gaze_status_codeRef;
    %if gaze is inside interactive canvas..
    elseif ( (cvsLeft_i < gaze_px(t,1)) && (gaze_px(t,1) < cvsRight_i) && (cvsTop_i < gaze_px(t,2)) && (gaze_px(t,2) < cvsBottom_i) )
      gaze_status(t,1) = cfg_gaze_status_codeInt;
    else
      gaze_status(t,1) = 0;
    endif
  endfor
endfunction

% returns string with jmol commands to select atoms of atom_index_array indexes, or return empty string if no index.
function current_jmol_script = jmol_scripting_selectAtom ( atom_index_array, transl_num)
  % if atom_index_array is not empty, fill it.
  if ~isempty(atom_index_array)
    current_jmol_script = "select ";
    % filling with the atom indexes
    for n=1:length(atom_index_array)
      current_jmol_script = strcat(current_jmol_script, " atomno=", num2str (atom_index_array(n)) );
      % write 'or' if this was not the last, or ';' if it was.
      if n < length(atom_index_array)
        current_jmol_script = strcat(current_jmol_script, " or");
      else
        current_jmol_script = strcat(current_jmol_script, ";");
      endif
    endfor
    current_jmol_script = strcat(current_jmol_script, [" color atoms TRANSLUCENT ", num2str(transl_num),";"] );
  # Else, leave it with empty space
  else
    current_jmol_script = "";
  endif
endfunction

% writes cell matrix with jmol commands for atoms transparency animation from heatmap_mw_int (or _ref). Needs to horzcat() to replay_rot_jmol_script
function replay_transp_jmol_script = replay_transparency (frame_count, atom_count, heatmap_mw)

  % building heatmap scale for animation
  for t=1:frame_count
    heatmap_mw_max(t) = max(heatmap_mw(t,:));
    heatmap_mw_min(t) = min(heatmap_mw(t,:));
    heatmap_mw_range(t) = heatmap_mw_max(t) - heatmap_mw_min(t);
    % compute scale from 1 (max) to nearly 0 for all atoms in each time frame
    for a=1:atom_count
      % 0.05 is minimal value to avoid division by 0
      heatscale_mw(t,a) = heatmap_mw(t,a) / max ( heatmap_mw_max(t), 0.05 );
    endfor
  endfor
  % build atom x frame matrix for "painting atoms with transparency" (scale 'e' of eights in jmol)
  transparency_matrix = 1 - round( heatscale_mw * 8)/8;
  replay_transp_jmol_script = cell(frame_count,9);
  % first frame, all atoms are set. loop through each transparency scale
  for e = 0:8
    atom_index_array = find ( transparency_matrix(1,:) == e/8 );
    current_jmol_script = jmol_scripting_selectAtom ( atom_index_array, e/8);
    replay_transp_jmol_script{1,e+1} = current_jmol_script;
  endfor
  % for the rest of the frames, and each scale
  for t = 2:frame_count
    for e = 0:8
      atom_index_array = find ( (transparency_matrix(t,:) == e/8 ) .* ( transparency_matrix(t-1,:) != transparency_matrix(t,:) ) );
      current_jmol_script = jmol_scripting_selectAtom ( atom_index_array, e/8);
      replay_transp_jmol_script{t,e+1} = current_jmol_script;
    endfor
  endfor

endfunction

% write file with jmol commands to animate the replay with gaze heatmap in 3D
function writeOutput_heatmapMw (filename, data_matrix_int, data_matrix_ref)
  tic(); printf("Writing heatmap output file..");
  %open file pointer
  xls_heatmapMw = xlsopen (filename, true);

  [xls_heatmapMw] = oct2xls (data_matrix_int, xls_heatmapMw, "jmol gaze int");
  [xls_heatmapMw] = oct2xls (data_matrix_ref, xls_heatmapMw, "jmol gaze ref");

  %close file pointer
  [xls_heatmapMw] = xlsclose (xls_heatmapMw);
  printf(".done!"); toc();
endfunction

%==========================
%SCRIPTS
% data checks for the slowest functions!
if and ( exist('raw_eyeT_data', 'var') == 0 , cfg_eyeT_input == false)
  warning("The eyeTracking data source is missing! Changing cfg_eyeT_input to true \n");
  cfg_eyeT_input = true;
endif
if and ( exist('raw_iRT_data', 'var') == 0 , cfg_iRT_input == false)
  warning("The iRT data source is missing! Changing cfg_iRT_input to true \n");
  cfg_iRT_input = true;
endif
% eyetracking data input (eyeT2oct)
if cfg_eyeT_input == true
  [raw_eyeT_data, raw_eyeT_header_data] = eyeT2oct (cfg_eyeT_input_filename);
else
  disp("Skipping eyeT file read.");
endif
% create eyeT_data and eyeT_header_data from respective raw data by adding epoch column to them
if cfg_add_epoch_from_milliseconds == true
  eyeT_data = horzcat( add_epoch (raw_eyeT_data, cfg_epoch_ms_col, cfg_epoch_anchor) , raw_eyeT_data );
  eyeT_header_data = horzcat ('EyeT Epoch', raw_eyeT_header_data);
else
  eyeT_data = raw_eyeT_data;
endif
% eyeTracking data pre-process: interpolation of missing data
if cfg_pupil_fix == true
  eyeT_data = interpolate_missing_data (eyeT_data, cfg_pupil_cols, cfg_pupil_missing_data);
endif
% iRT_data and session_data input (calculates session_row with either true or false)
if cfg_iRT_input == true
  [session_data,raw_iRT_data] = iRT2oct (cfg_iRT_input_filename);
  session_row = get_session_row (session_data, cfg_iRT_sessionID, cfg_iRT_taskID);
else
  disp("Skipping iRT file read.");
  session_row = get_session_row (session_data, cfg_iRT_sessionID, cfg_iRT_taskID);
endif
% iRT - eyeTracking data merge. Headers included
if cfg_data_merge == true
  % slice out desired iRT_data from raw_iRT_data (set desired columns in cfg_iRT_cols setting)
  [task_data, iRT_header_data] = slice_task_data (raw_iRT_data, cfg_iRT_sessionID, cfg_iRT_taskID, cfg_iRT_cols);

  % add column with angle distance from reference to plot resolugram
  frame_count = size (Q,1);
  Q_ref = cell2mat ( session_data(session_row,cfg_atom_matrix_ref_cols) );
  Q(:,1:4) = task_data(:,cfg_atom_matrix_quat_cols);
  resolugram = zeros ( frame_count, 1 );
  for t=1:frame_count
    resolugram(t) = abs( acos (Q(t,4)) - acos(Q_ref(4)) ) * ( 180/ pi() ) ;
  endfor
  iRT_header_data = horzcat(iRT_header_data, "resolugram dist");
  task_data = horzcat ( task_data, resolugram );

  % merge eyeT data into iRT data
  task_data = merge_data(eyeT_data, task_data, cfg_eyeT_cols);
  % merge headers
  merged_header_data = horzcat(iRT_header_data, eyeT_header_data(cfg_eyeT_cols));

  % plot resolugram
  if cfg_plot_resolugram == true
    figure (1);
    plot (0.1*(1:frame_count) , resolugram);
    title (strcat ("Resolugram") );
    axis ([ 0 frame_count*0.1 0 180 );
    xlabel("Task duration"); ylabel("Distance in degrees");
  endif
endif
% building interaction replay animation from temporal quaternion t x 4 array
replay_rot_jmol_script = cell(frame_count,1);
replay_rot_jmol_script{1} = strcat("moveto 0.0 QUATERNION {", num2str(Q(1,1:4)) , "};");
for t=2:frame_count
  if isequal( Q(t,:), Q(t-1,:))
    %if no rotation was made
    replay_rot_jmol_script{t} = "delay 0.1;";
  else
    %if some rotation was made
    replay_rot_jmol_script{t} = strcat("moveto 0.1 QUATERNION {", num2str(Q(t,1:4)) , "};");
  endif
endfor


% WRITE output file
if cfg_write_output == true
  writeOutput (cfg_output_filename, merged_header_data, task_data);
endif
% xyz data input
if cfg_xyz_input == true
  % get model_file_name used in the chosen session/task
  model_file_name = cell2mat( session_data(session_row,cfg_xyz_col) );
  % read file content
  [atom_count, atom_elem, atom_xyz] = get_xyz_data (strcat ("models/", model_file_name) );
  % normalize atom_xyz center of rotation to 0,0,0
  atom_xyz = normalize_jmol_rot_center (atom_xyz);
  %3D scatter of vertices/atoms without rotation
  if cfg_xyz_plot==true
    atom_cor = generate_color_vector (atom_count, atom_xyz, atom_elem);
    figure (2);
    scatter3 (atom_xyz(:,1), atom_xyz(:,2), atom_xyz(:,3), atom_cor(:,1), atom_cor(:,2:4));
    title ("3D Scatter of vertices");
    axis ("equal");
    xlabel("x"); ylabel("y"); zlabel("z");
  endif
endif


% compute atom matrices (temporal and reference) from xyz coordinates rotated acording to rotation data in chosen session
if cfg_atom_matrix == true
  % create temporal atom matrix of the interactive model
  atomInt_xyzRot = rotate_atom_xyz (task_data, cfg_atom_matrix_quat_cols, atom_xyz);
  % create atom matrix of the reference model
  atomRef_xyzRot = rotate_ref_atom_xyz (session_data, session_row, cfg_atom_matrix_ref_cols, atom_xyz);

  % compute xy pixel coordinate of canvas center (both reference and interactive)
  cvsInt_center = get_cvs_center (session_data, session_row, cfg_gaze_cvsInt_cols);
  cvsRef_center = get_cvs_center (session_data, session_row, cfg_gaze_cvsRef_cols);
  % get pixel to angstrom ratio of chosen row in session_data
  pxAngs_rate = session_data{session_row, cfg_gaze_pxAngs_rate_col};
  % get px coordinates of atoms xy projection (temporal/interactive and static/reference). (0,0) is top-left corner of screen
  atom_int_px = (atomInt_xyzRot(:,1:2,:) * pxAngs_rate) + cvsInt_center;
  atom_ref_px = (atomRef_xyzRot(:,1:2) * pxAngs_rate) + cvsRef_center;
  %plot above matrix in 2D at time set as cfg_atom_matrix_plot_t
  if and ( cfg_atom_matrix_plot == true, cfg_atom_matrix_plot_t > size(atomInt_xyzRot,3) )
    warning("The chosen frame index %i is outside the matrix range. Choose a reasonable frame index");
  elseif cfg_atom_matrix_plot == true
    atom_cor = generate_color_vector (atom_count, atom_xyz, atom_elem);
    figure (3);
    scatter (atomInt_xyzRot(:,1,cfg_atom_matrix_plot_t), atomInt_xyzRot(:,2,cfg_atom_matrix_plot_t), atom_cor(:,1), atom_cor(:,2:4));
    title ( strcat ("2D Scatter of vertices in interactive model at time=", num2str(cfg_atom_matrix_plot_t) ) );
    axis ("equal");
    xlabel("x"); ylabel("y");
  endif
  %plot above matrix in 2D
  if cfg_atom_matrix_ref_plot == true
    atom_cor = generate_color_vector (atom_count, atom_xyz, atom_elem);
    figure (4);
    scatter (atomRef_xyzRot(:,1), atomRef_xyzRot(:,2), atom_cor(:,1), atom_cor(:,2:4));
    title ("2D Scatter of vertices in reference model");
    axis ("equal");
    xlabel("x"); ylabel("y");
  endif

endif
% compute matrices of screen position distances between gazepoint and each atom
if cfg_gaze_dist_matrix == true
  gaze_px = task_data (:,cfg_gaze_cols);
  [gaze_atomInt_dist, gaze_atomRef_dist] = get_gaze_atom_dist (gaze_px, atom_int_px, atom_ref_px);
endif
% compute gaze_status in relation to canva : 2= int; 1= ref, 0= outside (TBD)
if cfg_gaze_status_array == true
  canvas_int = cell2mat (session_data (session_row,cfg_gaze_cvsInt_cols) ); %top, right, bottom, left side positions.
  canvas_ref = cell2mat (session_data (session_row,cfg_gaze_cvsRef_cols) ); %top, right, bottom, left side positions.
  %gaze_status = zeros ( size(gaze_px,1), 1); %n x 1
  gaze_status(:,1) = fill_gaze_status (gaze_px, canvas_ref, cfg_gaze_status_codeRef, canvas_int, cfg_gaze_status_codeInt);

endif
% Calculate heatmap from time spent in proximity of gaze during entire task
if cfg_gaze_heatmap == true
  printf("Calculating transparency gradient.");tic();

  %frame_count, atom_count, gaze_atomInt_dist
  frame_count = size(gaze_px,1);
  %initialize
  atom_gaze_alfa = ref_atom_gaze_alfa = zeros (atom_count,1);
  dist_ref = dist_int = zeros (frame_count,2,atom_count);
  dist_ref_exp = dist_int_exp = zeros (frame_count, atom_count);
  for a=1 : atom_count
%    for t=1 : size (Q,1)
    for t=1 : frame_count
      dist_int(t,1:2,a) = [ ( atom_int_px(a,1,t) - gaze_px(t,1) ).^2 , ( atom_int_px(a,2,t) - gaze_px(t,2) ).^2 ];
      dist_ref(t,1:2,a) = [ ( atom_ref_px(a,1) - gaze_px(t,1) ).^2 , ( atom_ref_px(a,2) - gaze_px(t,2) ).^2 ];
    endfor
  endfor
  for c=1:10
    for a=1 : atom_count
      cfg_gaussian_wdt = c*20;
      %gaussian formula: integral ( exp( - ( (x(t)-cx)^2 + (y(t)-cy)^2)/ (2*cfg_gaussian_wdt^2) ) dt)
      % compute and multiply by binary (0/1), testing if gaze is inside the canvas
      dist_ref_exp = exp (-(dist_ref(:,1,a)+dist_ref(:,2,a))/(2*cfg_gaussian_wdt^2) ) ;
      dist_int_exp = exp (-(dist_int(:,1,a)+dist_int(:,2,a))/(2*cfg_gaussian_wdt^2) ) ;
      %sum all
      ref_atom_gaze_alfa(a) = sum(dist_ref_exp.*(gaze_status==cfg_gaze_status_codeRef) );
      atom_gaze_alfa(a) = sum(dist_int_exp.*(gaze_status==cfg_gaze_status_codeInt) );
    endfor
    %plot(atom_gaze_alfa);
    xlswrite ("lista_alfas_atomos.xlsx", [cfg_gaussian_wdt;0;ref_atom_gaze_alfa], strcat("ref_",cfg_iRT_taskID), strcat (char (c+64), "2") );
    xlswrite ("lista_alfas_atomos.xlsx", [cfg_gaussian_wdt;0;atom_gaze_alfa], strcat(cfg_iRT_taskID), strcat (char (c+64), "2") );
%    printf("transparency gradient files written\n");
  endfor
  printf("concluido!");toc();
endif

% MOVING WINDOW HEATMAP

% Calculate a moving window heatmap from time spent in proximity of gaze during for the entire task
if cfg_gaze_heatmap_window == true
  printf("Calculating transparency gradient animation (writing may take some minutes):\n");tic();
  %initialize
  heatmap_mw_int = heatmap_mw_ref = zeros (atom_count,1); #(a,1)
  distMw_ref = distMw_int = zeros (frame_count,2,atom_count);  #(t,1:2,a)
%  distMw_ref_exp = distMw_int_exp = zeros (frame_count, atom_count); #(t,a)
  exp_distMw_ref = exp_distMw_int = zeros (frame_count, atom_count); #(t,a)
  ##gaussian formula: integral ( exp( - ( (x(t)-cx)^2 + (y(t)-cy)^2)/ (2*cfg_gaussian_wdt^2)) dt)

  for a=1 : atom_count
    for t=1 : frame_count
%      distMw_ref(t,1:2,a) = [ ( atom_ref_px(a,1)   - gaze_px(t,1) ).^2 , ( atom_ref_px(a,2)   - gaze_px(t,2) ).^2 ];
%      distMw_int(t,1:2,a) = [ ( atom_int_px(a,1,t) - gaze_px(t,1) ).^2 , ( atom_int_px(a,2,t) - gaze_px(t,2) ).^2 ];
      distMw_ref(t,a) = [ ( atom_ref_px(a,1)   - gaze_px(t,1) ).^2 + ( atom_ref_px(a,2)   - gaze_px(t,2) ).^2 ];
      distMw_int(t,a) = [ ( atom_int_px(a,1,t) - gaze_px(t,1) ).^2 + ( atom_int_px(a,2,t) - gaze_px(t,2) ).^2 ];
    endfor
  endfor

%  scatter (( atom_int_px(a,1,100) .- gaze_px(100,1) ),( atom_int_px(a,2,100) .- gaze_px(100,2) ) );

%  for c=1:10
  Q(:,1:4) = task_data(:,cfg_atom_matrix_quat_cols);
  c=1;
  for a=1 : atom_count
    % fill gaussian integral values
    exp_distMw_ref(:,a) = exp ( - distMw_ref(:,a) / (2*cfg_gaussian_wdt^2) ) ;
    exp_distMw_int(:,a) = exp ( - distMw_int(:,a) / (2*cfg_gaussian_wdt^2) ) ;
  endfor

  % for each atom, sums all values inside the moving window ("integrate") and register
  for t=1 : frame_count %building transparency gradient table with moving window
    %cfg_heatmap_mw_frame_length rule (spans from 's' to 't')
    s = max([t-cfg_heatmap_mw_frame_length, 1]);
    heatmap_mw_ref(t,1:atom_count) = sum (exp_distMw_ref(s:t,1:atom_count).*(gaze_status(s:t)==cfg_gaze_status_codeRef) );
    heatmap_mw_int(t,1:atom_count) = sum (exp_distMw_int(s:t,1:atom_count).*(gaze_status(s:t)==cfg_gaze_status_codeInt) );
  endfor



  replay_transp_jmol_script_int = horzcat(replay_rot_jmol_script, replay_transparency (frame_count, atom_count, heatmap_mw_int) );
  replay_transp_jmol_script_ref = replay_transparency (frame_count, atom_count, heatmap_mw_ref);
  printf("Writing jmol commands file (will take some time)...");tic();
  writeOutput_heatmapMw (cfg_heatmap_mw_filename, replay_transp_jmol_script_int, replay_transp_jmol_script_ref);
  printf("done!"); toc();


endif





%clear cfg variables for cleaner debugging
%clear cfg*
