pkg load io
%{
This script atempts to read the eyeTracker and iRT files to process them
and then merge both data in time. We used the UNIX EPOCH time to do this, as it
is the innate time used in most programming languages, especially javascript and
HTML. If the data from your eyetracker device/software does not contain this
kind of information, then it must be added using the epoch data produced from
iRT, or changing from another epoch system into the unix epoch.

Maybe usefull links:
https://wiki.octave.org/IO_package
%}

% clear all but the slowest variables to obtain
clear -exclusive *_data *_safe;

% SETTINGS
% set 'true' to let the script calculate (or recalculate/read/write again) the described values

% obtain needed values from iRT data. Mandatory
cfg_iRT_process = true; %SHOULD ALWAYS BE 'TRUE'
cfg_iRT_sessionID = 1682707472090; %session ID of the desired task
cfg_iRT_taskID = "bolaBastao_c"; %task ID of the desired task. Avoid using unsupported symbols for file names ( /\?|: )
%sample sessionIDs:    1682699553789    1682707472090    (sbj me)
%sample taskIDs:    bolaBastao_c    poligonFill    mrt

%{
plot_angDisp_multi (angDisp1_safe, "S1, task: bolaBastao", 15, [3,2,1])
plot_angDisp_multi (angDisp2_safe, "S2, task: bolaBastao", 15, [3,2,2])
plot_angDisp_multi (angDisp3_safe, "S1, task: poligonFill", 15, [3,2,3])
plot_angDisp_multi (angDisp4_safe, "S2, task: poligonFill", 15, [3,2,4])
plot_angDisp_multi (angDisp5_safe, "S1, task: mrt", 15, [3,2,5])
plot_angDisp_multi (angDisp6_safe, "S2, task: mrt", 15, [3,2,6])
%}

% read iRT .xlsx input file
cfg_iRT_input = true; %slow process. Set true to read a new file
cfg_iRT_input_filename = "iRT_data.xlsx"; %name of the iRT sheets file unpackaged by unpacking_sheets.m

% read eyeTracking .xlsx input file
cfg_eyeT_input = true; %slow process. Set true to read a new file
cfg_eyeT_input_filename = "raw_eyeT_1682707472090.xlsx"; %name of the input file (.xlsx, numbers only, no commas for decimals)

% fix missing pupil data by linear interpolation (best to always leave on with a new data array)
cfg_interpolate_missingVal = true;
cfg_interpolate_cols = [8,9]; % eyeT_data columns that need the interpolation. Pulil diameter, left/right, etc.
cfg_interpolate_vals = [0,-1]; % possible results of missing data to be identified and corrected

% compute and write file with table of jmol commands for the replay animation
cfg_replay_animation = true;
cfg_replay_animation_filename = ["output_copy_to_jmol_console ",num2str(cfg_iRT_sessionID)," ",cfg_iRT_taskID,".xlsx"];
cfg_plot_angDisp = true; % DRAW angular disparity plot (distance between target and interactive models, in degrees) (figure#1)

% merge eyeTracker data to the chosen iRT task (make sure the files are from the same session!)
cfg_data_merge = true;
cfg_iRT_cols = [3:8]; %range of desired data columns from raw_iRT_data. 5:8 is quaternion data, 3 is unix epoch
cfg_eyeT_cols = [1,2,4,6:9]; %range of desired data columns from eyeT_data. epoch data should be the 1st column
cfg_write_merge_output = true; % WRITE output file from task_data

% read .xyz file with atom data from the used model based on values in session_data
cfg_xyz_input = true;
cfg_xyz_col = 11; % index of the column to look for the modelName value in session_data
cfg_plot_xyz = false; % DRAW 3D vertices of the array of atoms colored acording to atom_elem (figure#2). Similar to what Jmol does, for debugging purposes.

% calculate temporal array of rotated atoms
cfg_atom_matrix = true;
cfg_atom_matrix_quat_cols = [3:6]; % column indexes of quaternions inside task_data array. Jmol quaternions order: [i j k r]
cfg_atom_matrix_tgt_cols = [7:10]; % column indexes for the target quaternion values inside session_data ("tgt_i","tgt_j","tgt_k","tgt_theta")
cfg_atom_matrix_tgt_plot = false; %2d scatter of the target model (figure#3)
cfg_atom_matrix_plot = false; %2d scatter of the interactive model at a set time (figure#4)
cfg_atom_matrix_plot_t = 1; %frame used in rotated cfg_atom_matrix_plot (figure#4)

% calculate gaze-canvas distances
cfg_gaze_dist_matrix = true;
cfg_gaze_cols = [11,12]; % column indexes of gaze x and y coordinates on screen. Should be in pixels, counting from top-left corner
cfg_gaze_scrSize_cols = [12,13]; % column indexes of screenSize values from session_data (width and height respectively)
cfg_gaze_cvsTgt_cols = [14,15,16,17]; % column indexes of target model canvas positions from session_data (in order: top, right, bottom, left)
cfg_gaze_cvsInt_cols = [18,19,20,21]; % column indexes of interactive model canvas positions from session_data (in order: top, right, bottom, left)
cfg_gaze_pxAngs_rate_col = [6]; %column index of pixels (screen distance) per angstrom (atomic distance unit in jmol) in session_data.

% calculate the status of the gaze in respect to where it is located: inside target canvas, interactive canvas, or neither
cfg_gaze_status_array = true;
cfg_gaze_status_codeTgt = 1; %condition in gaze_status, meaning that gaze was within Target model canvas
cfg_gaze_status_codeInt = 2; %condition in gaze_status, meaning that gaze was within Interactive model canvas
cfg_plot_angDisp_gaze_status = true; %plot the angular disparity with the line color based in the registered gaze_status

% calculate temporal transparency heatmap in 3D
cfg_gaze_heatmap_window = true;
cfg_heatmap_mw_frame_length = 20; %moving window length in frames for heatmap computation
cfg_gaussian_wdt = 50; %gaussian width in screen pixels, used in heatmap calculation

%{
   #=========================================#
   # DON'T MODIFY ANYTHING BELLOW THIS LINE! #
   #          (or do modify it if)           #
   #      (you know what you are doing)      #
   #=========================================#
%}

%FUNCTIONS
% recognize input file format and reads eyeTracker data (eyeT_data)
function [raw_eyeT_data, raw_eyeT_header_data] = eyeT2oct (filename)
  %open file pointer
  printf("Eyetracker data input: opening.."); tic();
  %raw_eyeT_data = xlsread (filename);
  xls_eyeT = xlsopen(filename);

  printf(".reading (lengthy step).");
  cell_data = xls2oct(xls_eyeT);
  printf(".parsing.");
  raw_eyeT_header_data = cell(cell_data(1,:));
  raw_eyeT_data = parsecell(cell_data);

  %closing file pointer
  printf(".closing.");
  [xls_eyeT] = xlsclose(xls_eyeT);
  printf(".OK! "); toc();
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
    printf("[1st valid row = %i].",first_valid_line);
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
  printf(".OK!\n");
endfunction

% read data from iRT .xlsx; 'session_arr' has session details, 'raw_arr' has the data bulk in a cell matrix
function [session_arr, raw_arr] = iRT2oct (filename)
  tic();
  printf("Reading iRT data (lengthy step)..");
  xls_iRT = xlsopen(filename);
  session_arr = xls2oct(xls_iRT,"sessions");
  raw_arr = xls2oct(xls_iRT,"data");
  [xls_iRT] = xlsclose(xls_iRT);
  printf(".OK!");
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
      printf(".[1st slice row: %i].", first_line);
      break;
    endif
  endfor
  % find last line of the slice
  last_line= -1;
  for n = first_line : size(raw_arr,1)
    %if n is from another task/session..
    if not ( and ( strcmp ( raw_arr{n,1}, num2str (session_ID) ), strcmp ( raw_arr{n,2}, task_ID) ) )
      last_line = n - 1;
      printf(".[last slice row: %i].", last_line);
      break;
    %if n is at the last row in file..
    elseif n == size(raw_arr,1)
      last_line = n;
      printf(".[last slice row: %i].", last_line);
      break;
    endif

  endfor
  % get slice of data
  numeric_data = cell2mat( raw_arr( first_line:last_line, desired_iRT_columns) );
  header_data = cell(raw_arr(1,desired_iRT_columns));
  printf(".OK! \n");
endfunction
% add column with angle disparity to plot angular disparity
function angDisp = compute_angDisp (Q, Q_tgt)
  angDisp = zeros ( size(Q,1), 1 );
  Q_tgt = Q_tgt / norm(Q_tgt);
  a_tgt = Q_tgt(1);
  b_tgt = Q_tgt(2);
  c_tgt = Q_tgt(3);
  d_tgt = Q_tgt(4);
  tgt_matrix_transposed = [d_tgt, c_tgt, -b_tgt, a_tgt; -c_tgt, d_tgt, a_tgt, b_tgt; b_tgt, -a_tgt, d_tgt, c_tgt; -a_tgt, -b_tgt, -c_tgt, d_tgt];
  % n x 4
  for t=1:size(Q,1)
    q(t,:) = Q(t,:) / norm ( Q(t,:) );
    r(t,:) = q(t,:) * tgt_matrix_transposed ;
    % q(4) should always be positive
    if r(t,4) < 0
      r(t,:) = -r(t,:);
    endif
    % compute angle disparity in radians
    angDisp(t) = 2 * acos(r(t,4));
    % transform radians to degrees
    angDisp(t) = angDisp(t) * 180 / pi;
  endfor
endfunction

% function to plot angular disparity
function plot_angDisp (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID)
  frame_count = size(angDisp,1);
  x = 0.1*(0:frame_count-1);
  y = angDisp;
  plot_angDisp_title = ["Resolugram - ",num2str(cfg_iRT_sessionID)," ",cfg_iRT_taskID];
  figure (1);
    plot (x , y);
    title (plot_angDisp_title);
    axis ([ 0 frame_count*0.1 0 180 ]);
    xlabel("Task duration"); ylabel("Distance in degrees");
endfunction

% function to plot angular diaparity with line colors (needs gaze_status values)
function plot_angDisp_colored_line (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, gaze_status)
  % computations
  frame_count = size(angDisp,1);
  plot_angDisp_title = ["Resolugram - ",num2str(cfg_iRT_sessionID)," ",cfg_iRT_taskID];
  x = 0.1*(0:frame_count-1);
  y = angDisp;
  color_0 = '#808080'; % gray
  color_1 = 'blue';
  color_2 = 'red';
  line_colors = {color_0, color_1, color_2};
  figure (10);
  hold on;
  % samples for the graph legend
  plot (x(1),y(1),'color', line_colors{1}, 'linewidth', 1.0 );
  plot (x(1),y(1),'color', line_colors{2}, 'linewidth', 1.0 );
  plot (x(1),y(1),'color', line_colors{3}, 'linewidth', 1.0 );

  % plot line with colored segments
  current_status = gaze_status(1);
  ln_len = 0;
  for i=2:frame_count
    ln_len = ln_len+1;
    if or ( i == frame_count, gaze_status(i) != current_status)
      %plot colored line segment
      plot ( x(i-ln_len:i), y(i-ln_len:i), 'color', line_colors{current_status+1}, 'linewidth', 1.0 );
      %prepare for next line segment
      current_status = gaze_status(i);
      ln_len = 0;
    endif
  endfor

  hold off;
  % Descriptions
  legend ('Outside', 'Target model', 'Interactive model');
  title (plot_angDisp_title);
  axis ([ 0 frame_count*0.1 0 180 ]);
  xlabel("Task duration"); ylabel("Distance in degrees");
  set(gca, 'ytick', 0:30:180);
endfunction


% function to plot angular disparity with background colors (needs gaze_status values)
function plot_angDisp_colored_bg (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, gaze_status)
  % computations
  frame_count = size(angDisp,1);
  plot_angDisp_title = ["Resolugram - ",num2str(cfg_iRT_sessionID)," ",cfg_iRT_taskID];
  x = 0.1*(0:frame_count-1);
  y = angDisp;
  color_0 = '#f0f0f0'; %outside
  color_1 = '#9090ff'; %target
  color_2 = '#ff9090'; %Interactive
  bg_colors = {color_0, color_1, color_2};
  figure (8);
  hold on;
  % samples for the graph legend
  plot (x(1),y(1),'color', bg_colors{1}, 'linewidth', 5.0 );
  plot (x(1),y(1),'color', bg_colors{2}, 'linewidth', 5.0 );
  plot (x(1),y(1),'color', bg_colors{3}, 'linewidth', 5.0 );

  % plot gray line with colored background
  current_status = gaze_status(1);
  ln_len = 0; % length of the line segment
  x_left = 0; % left corner of each rectangle
  for i=2:frame_count
    ln_len = ln_len+1;
    if i == frame_count || gaze_status(i) != current_status
      %plot colored line segment
      rectangle ( 'Position', [ (i-1-ln_len)/10 , 0, ln_len/10, 180], 'FaceColor', bg_colors{current_status+1}, 'EdgeColor', 'none');
      %prepare for next line segment
      current_status = gaze_status(i);
      ln_len = 0;
    endif
  endfor
  % plot the angular distance line
  plot ( x, y, 'color', '#000000', 'linewidth', 1.0 );
  hold off;
  % Descriptions
  legend ('Outside', 'Target model', 'Interactive model');
  title (plot_angDisp_title);
  axis ([ 0 frame_count*0.1 0 180 ]);
  xlabel("Task duration"); ylabel("Angular Disparity");
  set(gca, 'ytick', 0:30:180);
endfunction

% function to plot multiple angular disparity graphs in grid
function plot_angDisp_multi (angDisp, plot_angDisp_title, fig_n, sub_p)
  frame_count = size(angDisp,1);
  x_duration = 0.1*(1:frame_count);

  % fig_n is a unique pointer to a plot. If two plots have the same pointer, one of them will be overwritten
  figure (fig_n);
    % subplot: rows, columns, index of selected plot
    subplot (sub_p(1),sub_p(2),sub_p(3));
    ax = plot (x_duration , angDisp);
    title(plot_angDisp_title);
    xlabel("Task duration");
    ylabel ("Angular Disparity");
    axis ([0,Inf, 0,180] );
    set(gca, 'ytick', 0:45:180);
endfunction

% function to plot multiple angular disparity graphs with more data
function plot_angDisp_multi_yy (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, extra_series, fig_n, sub_p)
  frame_count = size(angDisp,1);
  x_duration = 0.1*(1:frame_count);
  plot_angDisp_title = ["Resolugram - ",num2str(cfg_iRT_sessionID)," ",cfg_iRT_taskID," + pupil data"];

  figure (fig_n);
    % subplot: rows, columns, index of selected plot
    subplot (sub_p(1),sub_p(2),sub_p(3));
    ax = plotyy (x_duration , angDisp, x_duration, extra_series);
    title(plot_angDisp_title);
    xlabel("Task duration");
    ylabel (ax(1), "Resolugram - Distance in degrees");
    ylabel (ax(2), 'Pupil Diameter');
    axis (ax(1), [0,Inf, 0,180] );
    axis (ax(2), 'autoy' ); % auto specifies the y-axis length
endfunction

% merge eyeT_data into iRT_data based on the nearest time values by nearest neighbours method
function mergedMatrix = merge_data (iRT_data, eyeT_data, desired_eyeT_columns)
  %Find nearest indexes in eyeT_data for each time point in iRT_data
  printf("Calculating nearest indexes.");
  eyeT_epoch_col = eyeT_data(:, 1); %eyeT epoch
  iRT_epoch_col = iRT_data(:, 1); %iRT epoch
  [~, nearest_indices] = min(abs(eyeT_epoch_col  - iRT_epoch_col'));
  % Merge desired eyeT_data columns into iRT_data
  printf(".merging to iRT_data.");
  mergedMatrix = [iRT_data, eyeT_data(nearest_indices, desired_eyeT_columns)];
  printf(".OK!\n");
endfunction

% write output file from merged data
function writeOutput_merged (filename, header, task_data, session_data, session_row)
  printf("Writing output file with merged data points.."); tic();
  % open file pointer
  xls_merged = xlsopen (filename, true);

  %write files
  merged_cell_matrix = vertcat (header, num2cell (task_data));
  [xls_merged] = oct2xls (merged_cell_matrix , xls_merged, "data");
  [xls_merged] = oct2xls ( vertcat (session_data(1,:), session_data(session_row,:) ) , xls_merged, "session");

  % close file pointer
  [xls_merged] = xlsclose ( xls_merged);
  printf(".OK!"); toc();
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
  printf(" OK! \n");
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

% rotate atom_xyz atom matrix based on the target model quaternion from session_data
function atomTgt_xyzRot = rotate_tgt_atom_xyz (session_data, session_row, cfg_atom_matrix_tgt_cols, atom_xyz)
  Q_tgt = cell2mat ( session_data(session_row,cfg_atom_matrix_tgt_cols) );
  for a = 1:size(atom_xyz,1);
    atomTgt_xyzRot(a,1:3) = ( rot_matrix (Q_tgt(1), Q_tgt(2), Q_tgt(3), Q_tgt(4)) * atom_xyz(a,1:3)' )' ;
  endfor
endfunction
% calculate, in pixels, the canvas center from the browser window
function canvas_center = get_cvs_center (session_data, row_index, cvs_pos_cols)
  canvas = cell2mat (session_data (row_index,cvs_pos_cols) ); %top, right, bottom, left side positions.
  canvas_center = [ (canvas(2) + canvas(4) )/2 , (canvas(1) + canvas(3) )/2 ];
endfunction
% compute pixel distance between gaze matrix and each atom of selected atom matrix
function [gaze_atomInt_dist, gaze_atomTgt_dist] = get_gaze_atom_dist (gaze_px, atom_int_px, atom_tgt_px)
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
        gaze_atomTgt_dist(t,a) = sqrt ( (atom_tgt_px(a,1)-gaze_px(t,1))^2 + (atom_tgt_px(a,2)-gaze_px(t,2))^2 );
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

% returns gaze_status (if gaze is inside canvas 1 or 2) from gaze pixel position and canvas extremities position in px [top right bottom left]
function gaze_status = fill_gaze_status (gaze_px, canvas_tgt, cfg_gaze_status_codeTgt, canvas_int, cfg_gaze_status_codeInt)
  frame_count = size(gaze_px, 1);
  %0,0 is top-left corner
  cvsTop_r = canvas_tgt(1);
  cvsRight_r = canvas_tgt(2);
  cvsBottom_r = canvas_tgt(3);
  cvsLeft_r = canvas_tgt(4);
  cvsTop_i = canvas_int(1);
  cvsRight_i = canvas_int(2);
  cvsBottom_i = canvas_int(3);
  cvsLeft_i = canvas_int(4);
  %fill gaze_status if gaze is within canvas bounduaries
  for t=1 : frame_count
    % if gaze is inside target canvas..
    if ( (cvsLeft_r < gaze_px(t,1)) && (gaze_px(t,1) < cvsRight_r) && (cvsTop_r < gaze_px(t,2)) && (gaze_px(t,2) < cvsBottom_r) )
      gaze_status(t,1) = cfg_gaze_status_codeTgt;
    %if gaze is inside interactive canvas..
    elseif ( (cvsLeft_i < gaze_px(t,1)) && (gaze_px(t,1) < cvsRight_i) && (cvsTop_i < gaze_px(t,2)) && (gaze_px(t,2) < cvsBottom_i) )
      gaze_status(t,1) = cfg_gaze_status_codeInt;
    else
      gaze_status(t,1) = 0;
    endif
  endfor
endfunction
% reads atom_index_array and returns jmol command string to change translucent property in indexed atoms. Returns empty if empty index.
function current_jmol_script = jmol_scripting_selectAtom ( atom_index_array, transl_num)
  % example of a non-empty return: "select atomno=2 or atomno=4 or atomno=29; color atoms TRANSLUCENT 0.375; "
  % if atom_index_array is not empty
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
    current_jmol_script = strcat(current_jmol_script, [" color atoms TRANSLUCENT ", num2str(transl_num, '%.3f'),";"] );
  # Else, leave it with empty space
  else
    current_jmol_script = "";
  endif
endfunction

% returns jmol commands cell matrix for atoms transparency animation from heatmap_mw_int (or _tgt). Needs to horzcat() to replay_int_jmol_script
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
function writeOutput_heatmapMw (filename, data_matrix_int, data_matrix_tgt)
  tic(); printf(["Writing jmol commands in file '",filename,"' for gaze heatmap animation replay .."]);
  %open file pointer
  xls_heatmapMw = xlsopen (filename, true);

  [xls_heatmapMw] = oct2xls (data_matrix_int, xls_heatmapMw, "jmol gaze int");
  [xls_heatmapMw] = oct2xls (data_matrix_tgt, xls_heatmapMw, "jmol gaze tgt");

  %close file pointer
  [xls_heatmapMw] = xlsclose (xls_heatmapMw);
  printf(".OK!"); toc();
endfunction

%==========================

%SCRIPTS
% data checks for the slowest functions!
if and ( exist('raw_iRT_data', 'var') == 0 , cfg_iRT_input == false)
  warning("The iRT data source is missing! Changed cfg_iRT_input to true for this execution.\n");
  cfg_iRT_input = true;
endif
if and ( exist('raw_eyeT_data', 'var') == 0 , cfg_eyeT_input == false)
  warning("The eyeTracking data source is missing! Changed cfg_eyeT_input to true for this execution.\n");
  cfg_eyeT_input = true;
endif
% eyetracking data input (eyeT2oct)
if cfg_eyeT_input == true
  [raw_eyeT_data, raw_eyeT_header_data] = eyeT2oct (cfg_eyeT_input_filename);
else
  disp("Skipping eyeT file read.");
endif
% eyeTracking data pre-process: interpolation of missing data
if cfg_interpolate_missingVal == true
  % error flag
  if any ( size ( raw_eyeT_data, 2 ) < cfg_interpolate_cols )
    error("Warning: The column selected is above the data size.");
  endif
  eyeT_data = interpolate_missing_data (raw_eyeT_data, cfg_interpolate_cols, cfg_interpolate_vals);
endif
% iRT_data and session_data input (calculates session_row, Q_tgt, Q and frame_count with either true or false)
if cfg_iRT_input == true
  [session_data,raw_iRT_data] = iRT2oct (cfg_iRT_input_filename);
else
  disp("Skipping iRT file read.");
endif

% obtain needed values from specified sessionID (cfg_iRT_sessionID) and taskID (cfg_iRT_taskID) inside session_data
if cfg_iRT_process == true
  % obtain row index of chosen task/session inside session_data table
  session_row = get_session_row (session_data, cfg_iRT_sessionID, cfg_iRT_taskID);
  % slice out desired iRT_data from raw_iRT_data (set desired columns in cfg_iRT_cols setting)
  [iRT_data, iRT_header_data] = slice_task_data (raw_iRT_data, cfg_iRT_sessionID, cfg_iRT_taskID, cfg_iRT_cols);

  % obtain quaternion coordinates array (1 x 4) of target model
  Q_tgt = cell2mat ( session_data(session_row,cfg_atom_matrix_tgt_cols) );
  % obtain quaternion coordinates matrix (t x 4) of interactive model
  Q(:,1:4) = iRT_data(:,cfg_atom_matrix_quat_cols);
  % obtain amount of rows/frames/data points from the duration of the executed task
  frame_count = size (Q,1);

  % compute angular disparity data
  angDisp = compute_angDisp (Q, Q_tgt);
  iRT_data = horzcat ( iRT_data, angDisp );
  if cfg_plot_angDisp == true
    plot_angDisp (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID);
  endif
endif

% iRT - eyeTracking data merge. Headers included
if cfg_data_merge == true
  % merge eyeT data into iRT data
  task_data = merge_data(iRT_data, eyeT_data, cfg_eyeT_cols);
  % merge headers
  merged_header_data = horzcat(iRT_header_data, "angDisp", raw_eyeT_header_data(cfg_eyeT_cols));

  % WRITE output file
  if cfg_write_merge_output == true
    writeOutput_merged (["outputMerge ", num2str(cfg_iRT_sessionID), " ", cfg_iRT_taskID, ".xlsx"], merged_header_data, task_data, session_data, session_row);
  endif

endif

% building interaction replay animation from temporal quaternion (t x 4 array)
if cfg_replay_animation == true
  %declare cell matrix
  replay_int_jmol_script = cell(frame_count,1);
  %fill first row (all opaque)
  replay_int_jmol_script{1} = strcat("moveto 0.0 QUATERNION {", num2str(Q(1,1:4)) , "};");
  % fill rest of rows
  for t=2:frame_count
    if isequal( Q(t,:), Q(t-1,:))
      %if no rotation was made, just delay 0.1 seconds
      replay_int_jmol_script{t} = "delay 0.1;";
    else
      %if some rotation was made, include it in the Jmol commands
      replay_int_jmol_script{t} = strcat("moveto 0.1 QUATERNION {", num2str(Q(t,1:4)) , "};");
    endif
  endfor
  tic(); printf("Writing file with jmol commands for rotation animation replay ..");
  %open file pointer
  xls_replay = xlsopen (cfg_replay_animation_filename, true);

  [xls_replay] = oct2xls (replay_int_jmol_script, xls_replay, "rotaion");

  %close file pointer
  [xls_replay] = xlsclose (xls_replay);
  printf(".OK!"); toc();
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
  if cfg_plot_xyz==true
    atom_cor = generate_color_vector (atom_count, atom_xyz, atom_elem);
    figure (2);
    scatter3 (atom_xyz(:,1), atom_xyz(:,2), atom_xyz(:,3), atom_cor(:,1), atom_cor(:,2:4));
    title ("3D Scatter of vertices");
    axis ("equal");
    xlabel("x"); ylabel("y"); zlabel("z");
  endif
endif


% compute atom matrices (temporal and target) from xyz coordinates rotated acording to rotation data in chosen session
if cfg_atom_matrix == true
  % create temporal atom matrix of the interactive model
  atomInt_xyzRot = rotate_atom_xyz (task_data, cfg_atom_matrix_quat_cols, atom_xyz);
  % create atom matrix of the target model
  atomTgt_xyzRot = rotate_tgt_atom_xyz (session_data, session_row, cfg_atom_matrix_tgt_cols, atom_xyz);

  % compute xy pixel coordinate of canvas center (both target and interactive)
  cvsInt_center = get_cvs_center (session_data, session_row, cfg_gaze_cvsInt_cols);
  cvsTgt_center = get_cvs_center (session_data, session_row, cfg_gaze_cvsTgt_cols);
  % get pixel to angstrom ratio of chosen row in session_data
  pxAngs_rate = session_data{session_row, cfg_gaze_pxAngs_rate_col};
  % get px coordinates of atoms xy projection (temporal/interactive and static/target). (0,0) is top-left corner of screen
  atom_int_px = (atomInt_xyzRot(:,1:2,:) * pxAngs_rate) + cvsInt_center;
  atom_tgt_px = (atomTgt_xyzRot(:,1:2) * pxAngs_rate) + cvsTgt_center;
  %plot interactive model atoms/vertices in 2D at given time set as cfg_atom_matrix_plot_t
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
  %plot target model atoms/vertices in 2D
  if cfg_atom_matrix_tgt_plot == true
    atom_cor = generate_color_vector (atom_count, atom_xyz, atom_elem);
    figure (4);
    scatter (atomTgt_xyzRot(:,1), atomTgt_xyzRot(:,2), atom_cor(:,1), atom_cor(:,2:4));
    title ("2D Scatter of vertices in target model");
    axis ("equal");
    xlabel("x"); ylabel("y");
  endif

endif
% compute matrices of screen position distances between gazepoint and each atom
if cfg_gaze_dist_matrix == true
  gaze_px = task_data (:,cfg_gaze_cols);
  [gaze_atomInt_dist, gaze_atomTgt_dist] = get_gaze_atom_dist (gaze_px, atom_int_px, atom_tgt_px);
endif
% compute gaze_status in relation to canva :  1= tgt; 2= int; 0= outside
if cfg_gaze_status_array == true
  canvas_int = cell2mat (session_data (session_row,cfg_gaze_cvsInt_cols) ); % top, right, bottom, left side positions.
  canvas_tgt = cell2mat (session_data (session_row,cfg_gaze_cvsTgt_cols) ); % top, right, bottom, left side positions.
  gaze_status = zeros ( size(gaze_px,1), 1); % n x 1
  gaze_status(:,1) = fill_gaze_status (gaze_px, canvas_tgt, cfg_gaze_status_codeTgt, canvas_int, cfg_gaze_status_codeInt);
endif
% plot angular disparity colored with gaze_status
if cfg_plot_angDisp_gaze_status == true
  plot_angDisp_colored_bg (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, gaze_status);
endif

% Calculate a moving window heatmap from atoms proximity of gaze in time during the entire task, (TBD: ponderated with a decay in time)
if cfg_gaze_heatmap_window == true
  printf("Calculating transparency gradient for replay animation (may take a while):\n");tic();
  %declaring variables
  heatmap_mw_int = heatmap_mw_tgt = zeros (atom_count,1); #(a,1)
  distMw_tgt = distMw_int = zeros (frame_count,2,atom_count);  #(t,1:2,a)
%  distMw_tgt_exp = distMw_int_exp = zeros (frame_count, atom_count); #(t,a)
  exp_distMw_tgt = exp_distMw_int = zeros (frame_count, atom_count); #(t,a)
  %gaussian formula: integral ( exp( - ( (x(t)-cx)^2 + (y(t)-cy)^2)/ (2*cfg_gaussian_wdt^2)) dt)

  for a=1 : atom_count
    for t=1 : frame_count
%      distMw_tgt(t,1:2,a) = [ ( atom_tgt_px(a,1)   - gaze_px(t,1) ).^2 , ( atom_tgt_px(a,2)   - gaze_px(t,2) ).^2 ];
%      distMw_int(t,1:2,a) = [ ( atom_int_px(a,1,t) - gaze_px(t,1) ).^2 , ( atom_int_px(a,2,t) - gaze_px(t,2) ).^2 ];
      distMw_tgt(t,a) = [ ( atom_tgt_px(a,1)   - gaze_px(t,1) ).^2 + ( atom_tgt_px(a,2)   - gaze_px(t,2) ).^2 ];
      distMw_int(t,a) = [ ( atom_int_px(a,1,t) - gaze_px(t,1) ).^2 + ( atom_int_px(a,2,t) - gaze_px(t,2) ).^2 ];
    endfor
  endfor

%  scatter (( atom_int_px(a,1,100) .- gaze_px(100,1) ),( atom_int_px(a,2,100) .- gaze_px(100,2) ) );

%  for c=1:10
  Q(:,1:4) = task_data(:,cfg_atom_matrix_quat_cols);
  c=1;
  for a=1 : atom_count
    % fill gaussian integral values
    exp_distMw_tgt(:,a) = exp ( - distMw_tgt(:,a) / (2*cfg_gaussian_wdt^2) ) ;
    exp_distMw_int(:,a) = exp ( - distMw_int(:,a) / (2*cfg_gaussian_wdt^2) ) ;
  endfor

  % for each atom, sums all values inside the moving window ("integrate") and register
  for t=1 : frame_count %building transparency gradient table with moving window
    %cfg_heatmap_mw_frame_length is the range in frames (0.1 second in 10Hz) for computing the relevance of each atom (spans from 's' to 't')
    s = max([t-cfg_heatmap_mw_frame_length, 1]);
    heatmap_mw_tgt(t,1:atom_count) = sum (exp_distMw_tgt(s:t,1:atom_count).*(gaze_status(s:t)==cfg_gaze_status_codeTgt) );
    heatmap_mw_int(t,1:atom_count) = sum (exp_distMw_int(s:t,1:atom_count).*(gaze_status(s:t)==cfg_gaze_status_codeInt) );
  endfor


  replay_tgt_jmol_script = repmat({"delay 0.1;"}, frame_count, 1);
  Q_tgt = cell2mat ( session_data(session_row,cfg_atom_matrix_tgt_cols) );
  replay_tgt_jmol_script{1} = ["moveto 0 QUATERNION {", num2str(Q_tgt(1:4)),"};"];
  replay_transp_jmol_script_tgt = horzcat( replay_tgt_jmol_script, replay_transparency (frame_count, atom_count, heatmap_mw_tgt) );
  replay_transp_jmol_script_int = horzcat( replay_int_jmol_script, replay_transparency (frame_count, atom_count, heatmap_mw_int) );

  writeOutput_heatmapMw (cfg_replay_animation_filename, replay_transp_jmol_script_int, replay_transp_jmol_script_tgt);
endif

