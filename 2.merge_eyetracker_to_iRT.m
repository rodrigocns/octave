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
cfg_eyeT_input = false; %slow, aprox. 1 min for 10k lines and 9 columns
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
cfg_iRT_input = false; %slow, aprox. 25 seconds for 4.5k lines and 8 columns
cfg_iRT_input_filename = "iRT_data.xlsx"; %name of the iRT sheets file unpackaged by unpacking_sheets.m

% merge eyeTracker data to the chosen iRT task (be sure that the files are from the same session)
cfg_data_merge = false;
cfg_iRT_sessionID = 1682707472090; %session ID of the desired task
cfg_iRT_taskID = "bolaBastao_c"; %task ID of the desired task
cfg_iRT_cols = [3:8]; %range of desired data columns from raw_iRT_data. 5:8 is quaternion data, 3 is unix epoch
cfg_eyeT_cols = [1,2,4,7:10]; %range of desired data columns from eyeT_data. epoch data was appended as the first column

% WRITE output file from task_data
cfg_write_output = false;
cfg_output_filename = strcat("mergeOutput_", cfg_iRT_taskID, num2str(cfg_iRT_sessionID), ".xlsx"); % name of the output file (.xlsx) from the merging of eyeT_data and iRT_data

% read .xyz file with atom data from the used model based on values in session_data
cfg_xyz_input = true;
cfg_xyz_col = 11; % index of the column to look for the modelName value in session_data
cfg_xyz_plot = false; % DRAW scatter3 of the array of atoms colored acording to atom_elem (image#1)

% calculate temporal array of rotated atoms
cfg_atom_matrix = true;
cfg_atom_matrix_quat_cols = [3:6]; % quaternion index of columns inside task_data array
cfg_atom_matrix_ref_cols = [7:10]; % column indexes for the reference quaternion values inside session_data ("ref_i","ref_j","ref_k","ref_theta")
cfg_atom_matrix_ref_plot = true; %2d scatter of the reference model (image#2)
cfg_atom_matrix_plot = true; %2d scatter of the interactive model at a set time (image#3)
cfg_atom_matrix_plot_t = 1; %frame used in rotated cfg_atom_matrix_plot (image#3)

% calculate gaze related things
cfg_gaze_matrix = false;
cfg_gaze_matrix_cols = [10,11]; % column indexes of gaze x and y coordinates on screen. Should be in pixels, counting from top-left corner
cfg_gaze_scrSize_cols = [12,13]; % column indexes of screenSize values from session_data (width and height respectively)
cfg_gaze_cvsRef_cols = [14,15,16,17]; % column indexes of reference model canvas positionsfrom session_data (in order: top, right, bottom, left)
cfg_gaze_cvsInt_cols = [18,19,20,21]; % column indexes of interactive model canvas positions from session_data (in order: top, right, bottom, left)
cfg_gaze_pxAngs_rate_col = [6]; %column index of pixels (screen distance) per angstrom (atomic distance unit in jmol) in session_data.

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
function atom_xyzRot = rotate_atom_xyz (task_data, atom_matrix_quat_cols, atom_xyz)
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
      atom_xyzRot(a,1:3,t) = (rot_vector(1:3,1:3,t)*atom_xyz(a,1:3)' )' ;
    endfor
  endfor
endfunction
% rotate atom_xyz atom matrix based on the reference model quaternion from session_data
function ref_atom_xyz = rotate_ref_atom_xyz (session_data, session_row, cfg_atom_matrix_ref_cols, atom_xyz)
  Q_ref = cell2mat ( session_data(session_row,cfg_atom_matrix_ref_cols) );
  for a = 1:size(atom_xyz,1);
    ref_atom_xyz(a,1:3) = ( rot_matrix (Q_ref(1), Q_ref(2), Q_ref(3), Q_ref(4)) * atom_xyz(a,1:3)' )' ;
  endfor
endfunction
% calculate, in pixels, the canvas center from the browser window
function canvas_center = get_cvs_center (session_data, row_index, cvs_pos_cols)
  canvas = cell2mat (session_data (row_index,cvs_pos_cols) ); %top, right, bottom, left side positions.
  canvas_center = [ (canvas(2) + canvas(4) )/2 , (canvas(1) + canvas(3) )/2 ];
endfunction
% compute temporal matrix of gaze distances to each atom
function [gazes] = gaze_calculations (task_data, gaze_matrix_cols, cvsRef_center, cvsInt_center, pxAngs_rate)

  %canvasRef = cell2mat (session_data (row_index,cvsRef_cols) );
  %canvasInt = cell2mat (session_data (row_index,cvsInt_cols) );
  %canvasRef_center = [ (canvasRef(2) + canvasRef(4) )/2 , (canvasRef(1) + canvasRef(3) )/2 ];
  %canvasInt_center = [ (canvasInt(2) + canvasInt(4) )/2 , (canvasInt(1) + canvasInt(3) )/2 ];

  %config_ref_center_px = [211,376]; %TBD
  %config_cvs_center_px = [615,379]; %TBD

  % read gaze xy position on screen in pixels (it NEEDS to be in pixels, starting at [0,0] in top-left )
  %gaze_px = task_data (:,gaze_matrix_cols);
  % compute relative position of gaze matrix as [0,0] being at the center of the respective canvas (int or ref)
  %gaze_ref_px = gaze_px - cvsRef_center;
  %gaze_int_px = gaze_px - cvsInt_center;
  %gaze_cvs_px = gaze_px - cvsInt_center;

  %atom_xy(:,1:2,:) = atom_xyzRot(:,1:2,:); %atom_xy(atoms, xy, frame) {centralized canvas}
  %atom_xy_px(:,:,:) = config_factor_px*atom_xy(:,:,:); %get px coordinates of atoms xy projection {already centralized}
  %ref_atom_xy_px(:,1:2) = ref_atom_xyz(:,1:2)*config_factor_px; %get pixel xy

  %get px coordinates of atoms xy projection centralized at the canvas (temporal/interactive and static/reference)
  %atom_xy_px = atom_xyzRot(:,1:2,:) * pxAngs_rate;
  %ref_atom_xy_px = ref_atom_xyz(:,1:2) * pxAngs_rate;


endfunction
% XXX
%==========================

%SCRIPTS
%data checks for the slowest functions!
if and ( exist('raw_eyeT_data', 'var') == 0 , cfg_eyeT_input == false)
  warning("The eyeTracking data source is missing! Changing cfg_eyeT_input to true \n");
  cfg_eyeT_input = true;
endif
if and ( exist('raw_iRT_data', 'var') == 0 , cfg_iRT_input == false)
  warning("The iRT data source is missing! Changing cfg_iRT_input to true \n");
  cfg_iRT_input = true;
endif
%eyetracking data input (eyeT2oct)
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
%eyeTracking data pre-process: interpolation of missing data
if cfg_pupil_fix == true
  eyeT_data = interpolate_missing_data (eyeT_data, cfg_pupil_cols, cfg_pupil_missing_data);
endif
%iRT data input (calculates session_row either way
if cfg_iRT_input == true
  [session_data,raw_iRT_data] = iRT2oct (cfg_iRT_input_filename);
  session_row = get_session_row (session_data, cfg_iRT_sessionID, cfg_iRT_taskID);
else
  disp("Skipping iRT file read.");
  session_row = get_session_row (session_data, cfg_iRT_sessionID, cfg_iRT_taskID);
endif
%iRT - eyeTracking data merge. Headers included
if cfg_data_merge == true
  % slice out desired iRT_data from raw_iRT_data (set desired columns in cfg_iRT_cols setting)
  [task_data, iRT_header_data] = slice_task_data (raw_iRT_data, cfg_iRT_sessionID, cfg_iRT_taskID, cfg_iRT_cols);
  % merge eyeT data into iRT data
  task_data = merge_data(eyeT_data, task_data, cfg_eyeT_cols);
  % merge headers
  merged_header_data = horzcat(iRT_header_data, eyeT_header_data(cfg_eyeT_cols));
endif
% WRITE output file
if cfg_write_output == true
  writeOutput (cfg_output_filename, merged_header_data, task_data);
endif
% xyz data input
if cfg_xyz_input == true
  % get model_file_name used in the chosen session/task
  session_row = get_session_row (session_data, cfg_iRT_sessionID, cfg_iRT_taskID);
  model_file_name = cell2mat( session_data(session_row,cfg_xyz_col) );
  % read file content
  [atom_count, atom_elem, atom_xyz] = get_xyz_data (strcat ("models/", model_file_name) );
  % normalize atom_xyz center of rotation to 0,0,0
  atom_xyz = normalize_jmol_rot_center (atom_xyz);
  %3D scatter of vertices/atoms without rotation
  if cfg_xyz_plot==true
    atom_cor = generate_color_vector (atom_count, atom_xyz, atom_elem);
    figure (1);
    scatter3 (atom_xyz(:,1), atom_xyz(:,2), atom_xyz(:,3), atom_cor(:,1), atom_cor(:,2:4));
    title ("3D Scatter of vertices");
    axis ("equal");
    xlabel("x"); ylabel("y"); zlabel("z");
  endif
endif


% create atom matrix (temporal and reference) from xyz coordinates rotated acording to quaternions in task_data
if cfg_atom_matrix == true
  atom_xyzRot = rotate_atom_xyz (task_data, cfg_atom_matrix_quat_cols, atom_xyz);
  %plot above matrix in 2D at time set as cfg_atom_matrix_plot_t
  if and ( cfg_atom_matrix_plot == true, cfg_atom_matrix_plot_t > size(atom_xyzRot,3) )
    warning("The chosen frame index %i is outside the matrix range. Choose a reasonable frame index");
  else
    atom_cor = generate_color_vector (atom_count, atom_xyz, atom_elem);
    figure (2);
    scatter (atom_xyzRot(:,1,cfg_atom_matrix_plot_t), atom_xyzRot(:,2,cfg_atom_matrix_plot_t), atom_cor(:,1), atom_cor(:,2:4));
    title ( strcat ("2D Scatter of vertices in interactive model at time=", num2str(cfg_atom_matrix_plot_t) ) );
    axis ("equal");
    xlabel("x"); ylabel("y");
  endif

  session_row = get_session_row (session_data, cfg_iRT_sessionID, cfg_iRT_taskID);
  % create atom matrix of the reference model
  ref_atom_xyz = rotate_ref_atom_xyz (session_data, session_row, cfg_atom_matrix_ref_cols, atom_xyz);

  %plot above matrix in 2D
  if cfg_atom_matrix_ref_plot == true
    atom_cor = generate_color_vector (atom_count, atom_xyz, atom_elem);
    figure (3);
    scatter (ref_atom_xyz(:,1), ref_atom_xyz(:,2), atom_cor(:,1), atom_cor(:,2:4));
    title ("2D Scatter of vertices in reference model");
    axis ("equal");
    xlabel("x"); ylabel("y");
  endif
endif

% gaze stuff calculation
if cfg_gaze_matrix == true
  % get row index of the chosen session
  session_row = get_session_row (session_data, cfg_iRT_sessionID, cfg_iRT_taskID);
  % compute xy pixel coordinate of canvas center (both reference and interactive)
  cvsRef_center = get_cvs_center (session_data, session_row, cfg_gaze_cvsRef_cols);
  cvsInt_center = get_cvs_center (session_data, session_row, cfg_gaze_cvsInt_cols);
  % read gaze xy position on screen in pixels (it NEEDS to be in pixels, starting at [0,0] in top-left )
  gaze_px = task_data (:,gaze_matrix_cols);
  % compute relative position of gaze matrix to the center of the respective canvas (int or ref)
  gaze_ref_px = gaze_px - cvsRef_center;
  gaze_int_px = gaze_px - cvsInt_center;

  % get pixel to angstrom ration of chosen session
  pxAngs_rate = session_data{session_row, cfg_gaze_pxAngs_rate_col};
  % get px coordinates of atoms xy projection centralized at the canvas (temporal/interactive and static/reference)
  atom_xy_px = atom_xyzRot(:,1:2,:) * pxAngs_rate;
  ref_atom_xy_px = ref_atom_xyz(:,1:2) * pxAngs_rate;

  %initialize arrays (time x atom_count)
  printf("Calculating gaze-atom distance array.."); tic();
  gaze_ref_atom_dist = gaze_atom_dist = zeros (size(Q,1),atom_count);
  for a=1:atom_count %for each atom
    for t=1 : size(Q,1) %for each point in time
      % compute distance
      %Formula: gaze_atom_dist = sqrt( (x-x')^2 + (y-y')^2 )
      gaze_atom_dist(t,a) = sqrt ( (atom_xy_px(a,1,t)-gaze_cvs_px(t,1))^2 + (atom_xy_px(a,2,t)-gaze_cvs_px(t,2))^2 );
      gaze_ref_atom_dist(t,a) = sqrt ( (ref_atom_xy_px(a,1)-gaze_cvs_px(t,1))^2 + (ref_atom_xy_px(a,2)-gaze_cvs_px(t,2))^2 );
    endfor
  endfor
  %time scan to fill in where are the gaze and its nearest atom
  for t=1 : size(Q,1)
    if ( -200<gaze_cvs_px(t,1) && gaze_cvs_px(t,1)<200 && -200<gaze_cvs_px(t,2) && gaze_cvs_px(t,2)<200  ) % marks in which part is the gaze
      gaze_status(t,1) = 2;
      [atom_closer_to_gaze(t,1),atom_closer_to_gaze(t,2)] = min (gaze_atom_dist(t,:));
    elseif ( -604<gaze_cvs_px(t,1) && gaze_cvs_px(t,1)<-204 && -203<gaze_cvs_px(t,2) && gaze_cvs_px(t,2)<197  )
      gaze_status(t,1) = 1;
      [atom_closer_to_gaze(t,1),atom_closer_to_gaze(t,2)] = min (gaze_ref_atom_dist(t,:));
    else
      gaze_status(t,1) = 0;
      atom_closer_to_gaze(t,1:2) = [0,0];
    endif
  endfor
  printf(".array calculated.");




endif


%clear cfg variables for cleaner debugging
%clear cfg*
