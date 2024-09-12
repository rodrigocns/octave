
pkg load io
%{
This script atempts to read the eyeTracker and iRT files to process them
and then merge both data in time. We used the UNIX EPOCH time to do this, as it
is the innate time used in most programming languages, especially javascript and
HTML. If the data from your eyetracker device/software does not contain this
kind of information, then it must be added using the epoch data produced from
iRT, or changing from another epoch system into the unix epoch.

For graphs, use Ctrl+F with "cfg_plot_"

Maybe usefull links:
https://wiki.octave.org/IO_package
%}

% clear all but the slowest variables to obtain
clear -exclusive *_data *_safe;

% SETTINGS
% set 'true' to let the script calculate (or recalculate/read/write again) the described values

% obtain needed values from iRT data. Mandatory
cfg_iRT_sessionID = "1682707472090"; %session ID of the desired task
cfg_iRT_taskID = "bolaBastao_c"; %task ID of the desired task. Avoid using unsupported symbols for file names ( /\?|: )
%sample sessionIDs:    1682699553789    1682707472090    (novice x expert)
%sample taskIDs:    bolaBastao_c    poligonFill    mrt

%{
plot_angDisp_multi (angDisp1_safe, "S1, task: bolaBastao", 15, [3,2,1])
plot_angDisp_multi (angDisp2_safe, "S2, task: bolaBastao", 15, [3,2,2])
plot_angDisp_multi (angDisp3_safe, "S1, task: poligonFill", 15, [3,2,3])
plot_angDisp_multi (angDisp4_safe, "S2, task: poligonFill", 15, [3,2,4])
plot_angDisp_multi (angDisp5_safe, "S1, task: mrt", 15, [3,2,5])
plot_angDisp_multi (angDisp6_safe, "S2, task: mrt", 15, [3,2,6])
%}

% 1.Temporal Data Input
% read iRT .xlsx input file
cfg_iRT_input = true; %slow process. Set true to read a new file
cfg_iRT_input_filename = "iRT data.xlsx"; %name of the iRT sheets file unpackaged by unpacking_sheets.m
cfg_iRT_filter = true; %SHOULD ALWAYS BE 'true'
cfg_plot_angDisp = false; % DRAW angular disparity plot (angular distance between target and interactive objects, in degrees) (figure#1)

% read eyeTracking .xlsx input file
cfg_eyeT_input = false; %slow process. Set true to read a new file
cfg_eyeT_input_filename = ["raw eyeT ", cfg_iRT_sessionID, ".xlsx"]; %name of the input file (.xlsx, numbers only, no commas for decimals)

% fix missing pupil data by linear interpolation (best to always leave on with a new data array)
cfg_interpolate_missingVal = true;
cfg_interpolate_cols = [8,9]; % eyeT_data columns that need the interpolation, such as Pulil diameter, etc.
cfg_interpolate_vals = [0,-1]; % possible values of missing data to be identified and corrected (if values match any of these, they will be corrected)

% 2.Data Synchronization
% merge eyeTracker data to the chosen iRT task (make sure the files are from the same session!)
cfg_data_merge = true;
cfg_iRT_cols = [3:8]; %range of desired data columns from raw_iRT_data. 5:8 is quaternion data, 3 is unix epoch
cfg_eyeT_cols = [1,2,4,6:9]; %range of desired data columns from eyeT_data. epoch data should be the 1st column
cfg_write_merge_output = false; % WRITE output file from task_data

% 3.Data Processing
%duration of a cycle (1/frequency)
cfg_t = 0.1;

% compute and write file with table of jmol commands for the replay animation
cfg_replay_animation = false; %ERROR?
cfg_replay_animation_filename = ["output jmol console ",num2str(cfg_iRT_sessionID)," ",cfg_iRT_taskID,".xlsx"];

% read .xyz file with atom data from the specified object in session_data
cfg_xyz_input = true;
cfg_xyz_col = 11; % index of the column to look for the modelName value in session_data
cfg_plot_xyz = false; % DRAW 3D vertices of the array of atoms colored acording to atom_elem (figure#2). Similar to what Jmol does, for debugging purposes.

% compute temporal array of rotated atoms
cfg_atom_matrix = true;
cfg_atom_matrix_quat_cols = [3:6]; % column indexes of quaternions inside task_data array. Jmol quaternions order: [i j k r]
cfg_atom_matrix_tgt_cols = [7:10]; % column indexes for the target quaternion values inside session_data ("tgt_i","tgt_j","tgt_k","tgt_theta")
cfg_atom_matrix_tgt_plot = false; %2d scatter of the target object (figure#3)
cfg_atom_matrix_plot = false; %2d scatter of the interactive object at a set time (figure#4)
cfg_atom_matrix_plot_t = 1; %frame used in rotated cfg_atom_matrix_plot (figure#4)

% compute gaze-canvas distances
cfg_gaze_dist_matrix = true;
cfg_gaze_cols = [11,12]; % column indexes of gaze x and y coordinates on screen. Should be in pixels, counting from top-left corner
cfg_gaze_scrSize_cols = [12,13]; % column indexes of screenSize values from session_data (width and height respectively)
cfg_gaze_cvsTgt_cols = [14,15,16,17]; % column indexes of target object canvas positions from session_data (in order: top, right, bottom, left)
cfg_gaze_cvsInt_cols = [18,19,20,21]; % column indexes of interactive object canvas positions from session_data (in order: top, right, bottom, left)
cfg_gaze_pxAngs_rate_col = [6]; %column index of pixels (screen distance) per angstrom (atomic distance unit in jmol) in session_data.
cfg_gaze_pxRatio_col = [22]; %column index of screen pixel ratio (the (atomic distance unit in jmol) in session_data.

% compute the region of the gaze in respect to where it is located: inside target canvas, interactive canvas, or neither
cfg_gaze_region_array = true;
cfg_gaze_region_codeTgt = 1; %condition in gaze_region, meaning that gaze was within Target object canvas
cfg_gaze_region_codeInt = 2; %condition in gaze_region, meaning that gaze was within Interactive object canvas
cfg_plot_angDisp_gaze_region = true; %plot angular disparity data with backgroud color based on the registered gaze_region
cfg_plot_gaze_region_mw = false; %plot angular disparity using moving average. -1 is tgt, 1 is int.
cfg_plot_angDisp_yy_pupil = false; %plot angular disparity data with pupil diameter data.

% compute temporal transparency heatmap in 3D
cfg_gaze_map = true;
cfg_gauss_wdt_screen = 50; % gaussian width (sigma) in screen pixels, used in heatmap calculation
cfg_gauss_wdt_time = 10; % gaussian width (sigma) for decay over time, in frames count (0.1s per frame, ).

% find the frame where exploration ends and adjustments begins
cfg_finding_stage_change = false;

% compute cumulative variation of angular disparity in time. Also discretize it.
cfg_angDisp_cumulative = true;
cfg_plot_angDisp_cumulative = false;

% compute the avg dist in time between center and recently gazed atoms (NO gazemap matrix)
cfg_gaze_dispersion = true;
cfg_plot_gaze_dispersion = false;
cfg_gaze_avg_dist = false;
cfg_gaze_window_size = 10;
cfg_gaze_centers_delta = true;  % compute distance between gaze centers
cfg_plot_gaze_centers_delta = false;
cfg_gaze_dispersion_center_speed = true; % Compute gaze center speed (deslocation / time)
cfg_plot_gaze_dispersion_center_speed = false;
%
cfg_gaze_movement = true;   % Compute gaze movement (linear, horizontal and vertical only)
cfg_plot_gaze_speed = false; % plot 4 graphs: saccade horizontal, vertical and absolute speeds (movement/time)+angDisp
cfg_plot_gaze_movement = false;    % plot x and y coordinates of gaze
cfg_plot_gaze_movement_cumulative = false;  %plot cumulative gaze movement
cfg_draw_atom_and_gaze = true; % plot both atom/vertices and gaze points in 2 graphs (target and interactive models)

%{
   #=========================================#
   # DON'T MODIFY ANYTHING BELLOW THIS LINE! #
   #       (  or do modify it if you )       #
   #       ( know what you are doing )       #
   #=========================================#
%}

%USER INPUT
%session ID
prpt_sessionID = inputdlg ("Insert session ID value. Ex: 1682699553789, 1682707472090", "Input Session ID");
if (isempty ( char (prpt_sessionID) ) == 1 )
  printf ( cstrcat ("Blank input, using default sessionID value: '",cfg_iRT_sessionID,"'.\n"));
else
  cfg_iRT_sessionID = char(prpt_sessionID);
  if cfg_iRT_sessionID == "novice" %redirect "novice" to its ID
    cfg_iRT_sessionID  = "1682699553789";
  endif
  printf ( cstrcat ("Using input value: '",cfg_iRT_sessionID,"'.\n"));
endif
%task ID
prpt_taskID = inputdlg ("Insert task ID value. Ex: bolaBastao\_c, poligonFill, mrt", "Input Task ID");
if (isempty ( char (prpt_taskID) ) == 1 )
  printf ( cstrcat ("Blank input, using default taskID value: '",cfg_iRT_taskID,"'.\n"));
else
  cfg_iRT_taskID = char(prpt_taskID);
  printf ( cstrcat ("Using input value: '",cfg_iRT_taskID ,"'.\n"));
endif

%iRT data filename
prpt_iRT_fname = inputdlg ("Insert interactive Rotations Task (iRT) unpackaged data filename. Ex: iRT data.xlsx", "Input iRT unpackaged data filename");
if (isempty ( char (prpt_iRT_fname) ) == 1 )
  printf ( cstrcat ("Blank input, using default iRT filename: '",cfg_iRT_input_filename,"'.\n"));
else
  cfg_iRT_input_filename = char(prpt_iRT_fname);
  printf ( cstrcat ("Using input value: '",cfg_iRT_input_filename ,"'.\n"));
endif

%eye tracker data filename
prpt_eyeT_fname = inputdlg ("Insert eye tracker data filename. Ex: raw eyeT 1682707472090.xlsx", "Input eye tracker data filename");
if (isempty ( char (prpt_eyeT_fname ) ) == 1 )
%  printf ( cstrcat ("Blank input, using default value: '",cfg_eyeT_input_filename,"'.\n"));
  cfg_eyeT_input_filename = ["raw eyeT ", char(cfg_iRT_sessionID), ".xlsx"];
  printf ( cstrcat ("Blank input, using default eye tracker filename: '",cfg_eyeT_input_filename,"'.\n"));
else
  cfg_eyeT_input_filename = char(prpt_eyeT_fname );
  printf ( cstrcat ("Using input value: '",cfg_eyeT_input_filename ,"'.\n"));
endif

%FUNCTIONS
% recognize input file format and reads eyeTracker data (eyeT_data)
function [raw_eyeT_data, raw_eyeT_header_data] = eyeT2oct (filename)
  %open file pointer
  printf("Eyetracker data input: opening.."); tic();
  %raw_eyeT_data = xlsread (filename); %commented out. Using other functions
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
      printf(".from row %i.", first_line);
      break;
    endif
  endfor
  % find last line of the slice
  last_line= -1;
  for n = first_line : size(raw_arr,1)
    %if n is from another task/session..
    if not ( and ( strcmp ( raw_arr{n,1}, num2str (session_ID) ), strcmp ( raw_arr{n,2}, task_ID) ) )
      last_line = n - 1;
      printf(".to %i].", last_line);
      break;
    %if n is at the last row in file..
    elseif n == size(raw_arr,1)
      last_line = n;
      printf(".to %i].", last_line);
      break;
    endif

  endfor
  % get slice of data
  numeric_data = cell2mat( raw_arr( first_line:last_line, desired_iRT_columns) );
  header_data = cell(raw_arr(1,desired_iRT_columns));
  printf(".OK! \n");
endfunction
% compute angular disparity column
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

% PLOT temporal data (such as angular disparity x time)
function plot_temporal_data (temporal_data, data_label, cfg_iRT_sessionID, cfg_iRT_taskID, axis_sizes=0, ytick_val=0, fig_num=0)
  % axis_sizes is last value of the data ([0,180] for a 0 to 180 degree)
  frame_count = size(temporal_data,1); %last x value
  x = 0.1*(0:frame_count-1); %x values (time)
  y = temporal_data;
%  temporal_data_title = [data_label," - ",num2str(cfg_iRT_sessionID)," ",cfg_iRT_taskID];
  temporal_data_title = [data_label];
  % format: fig_num = [lines, cols, subplot id, figure id]. Ex: [3,1,2,24]
  if sum(fig_num) == 0
    figure;
  else
    if length (fig_num) == 4
      figure ( fig_num(4) );
    endif
    subplot ( fig_num(1) , fig_num(2) , fig_num(3) );
  endif
  plot (x , y);
    title (temporal_data_title);
    if axis_sizes == 0 %if all axis_sizes values are 0 (default), ignore it.
      axis ([ 0, frame_count*0.1]);
    else %
      axis ([ 0, frame_count*0.1, axis_sizes]);
    endif
    xlabel("Task duration"); ylabel(data_label);
    if ytick_val != 0;
      set(gca, 'ytick', ytick_val);
      printf("ytick set as 'ytick_val'\n");%debug
    endif
endfunction
% PLOT "non-temporal" data (such as gaze ref x gaze int). TODO: update all plot_data_OLD calls to the new one
function plot_data_OLD (x_data,x_label, y_data,y_label, cfg_iRT_sessionID, cfg_iRT_taskID, axis_size=0, ytick_val=0, fig_num=0)
  % here, axis_size is a 1x2 or 1x4 array for the graph's axis property
  frame_count = size(y_data,1); %REMOVE
  x = x_data; %x values (NOT time)
  y = y_data;
  data_title = [y_label," x ",x_label," - ",num2str(cfg_iRT_sessionID)," ",cfg_iRT_taskID];
  % format: fig_num = [lines, cols, subplot id, figure id]. Ex: [3,1,2,24]
  if sum(fig_num) == 0
    figure;
  else
    if length (fig_num) == 4
      figure ( fig_num(4) );
    endif
    subplot ( fig_num(1) , fig_num(2) , fig_num(3) );
  endif
    plot (x , y);
    title (data_title);
    if axis_size != 0 %if axis_sizes values are not 0 (default), use them.
      axis (axis_size);
    endif
    xlabel(x_label); ylabel(y_label);
    if ytick_val != 0;
      set(gca, 'ytick', ytick_val);
      printf("ytick set as 'ytick_val'\n");%debug
    endif
endfunction
% PLOT "non-temporal" data (such as gaze ref x gaze int)
function plot_data (x_data,x_label, y_data,y_label, axis_size=0, ytick_val=0, fig_num=0)
  % here, axis_size is a 1x2 or 1x4 array for the graph's axis property
  frame_count = size(y_data,1); %REMOVE
  x = x_data; %x values (NOT time)
  y = y_data;
  data_title = [y_label," x ",x_label];
    % format: fig_num = [lines, cols, subplot id, figure id]. Ex: [3,1,2,24]
  if sum(fig_num) == 0
    figure;
  else
    if length (fig_num) == 4
      figure ( fig_num(4) );
    endif
    subplot ( fig_num(1) , fig_num(2) , fig_num(3) );
  endif
    plot (x , y);
    title (data_title);
    if axis_size != 0 %if axis_sizes values are not 0 (default), use them.
      axis (axis_size);
    endif
    xlabel(x_label); ylabel(y_label);
    if ytick_val != 0;
      set(gca, 'ytick', ytick_val);
      printf("ytick set as 'ytick_val'\n");%debug
    endif
endfunction
% PLOT angular disparity with background colors
function plot_angDisp_colored_bg (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, gaze_region, sub= [0,0])
  % decide if create new figure or insert inside subplot
  if sub(1) == 0
    %common plot if sub is 0 (default).
    figure;
  else
    % Else, inserts in existing subplot of sub(1) lines and single column, at position sub(2).
    subplot(sub(1),1,sub(2));
  endif

  % defining colors for bg
  color_0 = '#f0f0f0'; %outside
  color_1 = '#9090ff'; %target
  color_2 = '#ff9090'; %Interactive
  color_line = '#000000'; % Angular Disparity line
  bg_colors = {color_0, color_1, color_2};
  legend_labels = {'Outside', 'Target model', 'Interactive model', 'Angular disparity'};
  % samples for the graph legend
  hold on;
  plot (NaN,NaN,'color', bg_colors{1}, 'linewidth', 5.0, 'DisplayName', legend_labels{1} ); %dummy areas
  plot (NaN,NaN,'color', bg_colors{2}, 'linewidth', 5.0, 'DisplayName', legend_labels{2} );
  plot (NaN,NaN,'color', bg_colors{3}, 'linewidth', 5.0, 'DisplayName', legend_labels{3} );
  plot (NaN,NaN, 'color', color_line, 'linewidth', 1.0 ); % dummy line

  % data computation
  frame_count = size(angDisp,1);
  plot_angDisp_title = ["Angular disparity - ",num2str(cfg_iRT_sessionID)," ",cfg_iRT_taskID];
  x = 0.1*(0:frame_count-1);
  y = angDisp;

  % plot colored background
  current_region = gaze_region(1);
  ln_len = 0; % length of the line segment
  x_left = 0; % left corner of each rectangle
  for i=2:frame_count
    ln_len = ln_len+1;
    if i == frame_count || gaze_region(i) != current_region
      %plot colored bg rectangle
      rectangle ( 'Position', [ (i-1-ln_len)/10 , 0, ln_len/10, 180], 'FaceColor', bg_colors{current_region+1}, 'EdgeColor', 'none');
      %prepare for next bg rectangle
      current_region = gaze_region(i);
      ln_len = 0;
    endif
  endfor

  % plot the angular distance line
  plot ( x, y, 'color', color_line, 'linewidth', 1.0 );

  hold off;
  % Descriptions
  legend ('Outside', 'Target model', 'Interactive model', 'Angular disparity');
  title (plot_angDisp_title);
  axis ([ 0 frame_count*0.1 0 180 ]);  % xi xf yi yf
  xlabel("Task duration"); ylabel("Angular Disparity");
  set(gca, 'ytick', 0:30:180);
endfunction
% PLOT any temporal data with background colors
function plot_temporal_data_colored_bg (temporal_data, data_label, cfg_iRT_sessionID, cfg_iRT_taskID, gaze_region, fig_num= 0)
  % decide if new figure or subplot. format: fig_num = [lines, cols, subplot id, figure id*]. Ex: [3,1,2,24]
  % *: empty value if figure already exists.
  if sum(fig_num) == 0
    figure;
  else
    if length (fig_num) == 4
      figure ( fig_num(4) );
    endif
    subplot ( fig_num(1) , fig_num(2) , fig_num(3) );
  endif

  % defining colors for bg
  color_0 = '#f0f0f0'; %outside
  color_1 = '#9090ff'; %target
  color_2 = '#ff9090'; %Interactive
  color_line = '#000000'; % data line
  bg_colors = {color_0, color_1, color_2};
  legend_labels = {'Outside', 'Target model', 'Interactive model', data_label};
  % samples for the graph legend
  hold on;
  plot (NaN,NaN,'color', bg_colors{1}, 'linewidth', 5.0, 'DisplayName', legend_labels{1} ); %dummy areas
  plot (NaN,NaN,'color', bg_colors{2}, 'linewidth', 5.0, 'DisplayName', legend_labels{2} );
  plot (NaN,NaN,'color', bg_colors{3}, 'linewidth', 5.0, 'DisplayName', legend_labels{3} );
  plot (NaN,NaN, 'color', color_line, 'linewidth', 1.0 ); % dummy line

  % data computation
  frame_count = size(temporal_data,1);
  plot_title = [data_label," - ",num2str(cfg_iRT_sessionID)," ",cfg_iRT_taskID];
  x = 0.1*(0:frame_count-1);
  y = temporal_data;

  % plot colored background
  current_region = gaze_region(1);
  ln_len = 0; % length of the line segment
  x_left = 0; % left corner of each rectangle
  for i=2:frame_count
    ln_len = ln_len+1;
    if i == frame_count || gaze_region(i) != current_region
      %plot colored bg rectangle
      rectangle ( 'Position', [ (i-1-ln_len)/10 , 0, ln_len/10, 180], 'FaceColor', bg_colors{current_region+1}, 'EdgeColor', 'none');
      %prepare for next bg rectangle
      current_region = gaze_region(i);
      ln_len = 0;
    endif
  endfor

  % plot the angular distance line
  plot ( x, y, 'color', color_line, 'linewidth', 1.0 );

  hold off;
  % Descriptions
  legend ('Outside', 'Target model', 'Interactive model', 'Angular disparity');
  title (plot_title);
  axis ([ 0 frame_count*0.1 0 180 ]);  % xi xf yi yf
  xlabel("Task duration"); ylabel("Angular Disparity");
  set(gca, 'ytick', 0:30:180);
endfunction
% colored bg for gaze_region_mw. TBD.
function plot_gazeRegion_colored_bg_mw (gaze_region_mw, cfg_iRT_sessionID, cfg_iRT_taskID, sub= [0,0])
  % decide if create new figure or insert inside subplot
  if sub(1) == 0
    %common plot if sub is 0 (default).
    figure;
  else
    % Else, inserts in existing subplot of sub(1) lines and single column, at position sub(2).
    %example: sub = [3,2] goes to subplot(3,1,2)
    subplot(sub(1),1,sub(2));
  endif

  % defining colors for bg
  color_0 = '#f0f0f0'; %outside
  color_1 = '#9090ff'; %target
  color_2 = '#ff9090'; %Interactive
  color_line = '#000000'; % Angular Disparity line
  bg_colors = {color_0, color_1, color_2};
  legend_labels = {'Outside', 'Target model', 'Interactive model', 'Angular disparity'};
  % samples for the graph legend
  hold on;
  plot (NaN,NaN,'color', bg_colors{1}, 'linewidth', 5.0, 'DisplayName', legend_labels{1} ); %dummy areas
  plot (NaN,NaN,'color', bg_colors{2}, 'linewidth', 5.0, 'DisplayName', legend_labels{2} );
  plot (NaN,NaN,'color', bg_colors{3}, 'linewidth', 5.0, 'DisplayName', legend_labels{3} );
  plot (NaN,NaN, 'color', color_line, 'linewidth', 1.0 ); % dummy line

  % data computation
  frame_count = size(angDisp,1);
  plot_angDisp_title = ["Gaze region moving avg. - ",num2str(cfg_iRT_sessionID)," ",cfg_iRT_taskID];
  x = 0.1*(0:frame_count-1);
  y = angDisp;
  % array like gaze_region for coloring gaze_region_mw data
  gaze_ids = zeros(size(gaze_region_mw));
  for i=1:length(gaze_region_mw)
    if gaze_region_mw(i) > 0
      gaze_ids(i)= 2;
    elseif gaze_region_mw(i) < 0
      gaze_ids(i) = 1;
    else
      gaze_ids(i) = 0;
    endif
  endfor

  % plot colored background
  current_region = gaze_ids(1);
  ln_len = 0; % length of the line segment
  x_left = 0; % left corner of each rectangle
  for i=2:frame_count
    ln_len = ln_len+1;
    if i == frame_count || gaze_ids(i) != current_region
      %plot colored bg rectangle
      rectangle ( 'Position', [ (i-1-ln_len)/10 , -1, ln_len/10, 1], 'FaceColor', bg_colors{current_region+1}, 'EdgeColor', 'none');
      %prepare for next bg rectangle
      current_region = gaze_ids(i);
      ln_len = 0;
    endif
  endfor

  % plot the angular distance line
  plot ( x, y, 'color', color_line, 'linewidth', 1.0 );

  hold off;
  % Descriptions
  legend ('Outside', 'Target model', 'Interactive model', 'Gaze region');
  title (plot_angDisp_title);
  axis ([ 0 frame_count*0.1 -1 1 ]); % xi xf yi yf
  xlabel("Task duration"); ylabel("Gaze region moving average");
  set(gca, 'ytick', -1:0.5:1);
endfunction
% colored bg for gaze_region_mw. Clean background for manual painting TBD.
function plot_gazeRegion_colored_bg_mw_simple (gaze_region_mw, cfg_iRT_sessionID, cfg_iRT_taskID, sub= [0,0])
  % decide if create new figure or insert inside subplot
  if sub(1) == 0
    %common plot if sub is 0 (default).
    figure;
  else
    % Else, inserts in existing subplot of sub(1) lines and single column, at position sub(2).
    subplot(sub(1),1,sub(2));
  endif

  % defining colors for bg
  color_0 = '#f0f0f0'; %outside
  color_1 = '#9090ff'; %target
  color_2 = '#ff9090'; %Interactive
  color_line = '#000000'; % Plot line
  bg_colors = {color_0, color_1, color_2};
  legend_labels = {'Gaze at target model', 'Gaze at interactive model', 'Gaze region moving avg.'};
  % samples for the graph legend
  hold on;
  plot (NaN,NaN,'color', bg_colors{2}, 'linewidth', 5.0, 'DisplayName', legend_labels{2} );
  plot (NaN,NaN,'color', bg_colors{3}, 'linewidth', 5.0, 'DisplayName', legend_labels{3} );
  plot (NaN,NaN, 'color', color_line, 'linewidth', 1.0 ); % dummy line

  % data computation
  frame_count = size(gaze_region_mw,1);
  plot_title = ["Gaze region moving avg. - ",num2str(cfg_iRT_sessionID)," ",cfg_iRT_taskID];
  x = 0.1*(0:frame_count-1);
  y = gaze_region_mw;
  % split gaze_region_mw data
  gaze_ids_pos = gaze_ids_neg = zeros(size(gaze_region_mw));
  for i=1:frame_count
    if gaze_region_mw(i) < 0
      gaze_ids_neg(i) = gaze_region_mw(i);
    elseif gaze_region_mw(i) > 0
      gaze_ids_pos(i) = gaze_region_mw(i);
    endif
  endfor


  % plot the angular distance line
  plot ( x, y, 'color', color_line);
  hold on;
  h = area(0.1*(0:frame_count-1),[gaze_ids_neg,gaze_ids_pos] );
  set(h(1), 'FaceColor', '#9090ff');
  set(h(2), 'FaceColor', '#ff9090');
  hold off;
  % Descriptions
  legend ('Target model', 'Interactive model', 'Gaze region moving avg.');
  title (plot_title);
  axis ([ 0 frame_count*0.1 -1 1 ]); % xi xf yi yf
  xlabel("Task duration"); ylabel("Gaze region moving average");
  set(gca, 'ytick', -1:0.50:1);
endfunction

% PLOT multiple angular disparity graphs in grid
function plot_angDisp_multi (angDisp, plot_title, fig_n, sub_p)
  frame_count = size(angDisp,1);
  x_duration = 0.1*(1:frame_count);

  % fig_n is a unique pointer to a plot. If two plots have the same pointer, one of them will be overwritten
  figure (fig_n);
    % subplot: rows, columns, index of selected plot
    subplot (sub_p(1),sub_p(2),sub_p(3));
    ax = plot (x_duration , angDisp);
    title(plot_title);
    xlabel("Task duration");
    ylabel ("Angular Disparity");
    axis ([0,Inf, 0,180] );
    set(gca, 'ytick', 0:45:180);
endfunction

% PLOT temporal graphs with two Y axis
function plot_temporal_yy(temporal_data,temporal_name, extra_data,extraName, cfg_iRT_sessionID, cfg_iRT_taskID, fig_num=0)
  % format: fig_num = [lines, cols, subplot id, figure id]. Ex: [3,1,2,24]
  if sum(fig_num) == 0
    figure;
  else
    if length (fig_num) == 4
      figure ( fig_num(4) );
    endif
    subplot ( fig_num(1) , fig_num(2) , fig_num(3) );
  endif

  %time management
  frame_count = size(temporal_data,1);
  x_duration = 0.1*( (1:frame_count)-1 );

  temporal_data_title = [temporal_name," + ",extraName," - ",num2str(cfg_iRT_sessionID)," ",cfg_iRT_taskID];
  [ax,h1,h2] = plotyy (x_duration , temporal_data, x_duration, extra_data);
  title(temporal_data_title);
  xlabel("Task duration");
  ylabel (ax(1), temporal_name);
  ylabel (ax(2), extraName);
  axis (ax(1), 'autoy' );
  axis (ax(2), 'autoy' ); % auto specifies the y-axis length
endfunction
% PLOT angular disparity graphs with more data
function plot_angDisp_yy_other (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, extra_series, seriesName, sub=[0,0])
  if sub(1) == 0
    %common plot if sub is 0 (default).
    figure;
  else
    % Else, inserts in subplot of sub(1) lines and single column, at position sub(2).
    subplot(sub(1),1,sub(2));
  endif

  frame_count = size(angDisp,1);
  x_duration = 0.1*( (1:frame_count)-1 );
  plot_angDisp_title = ["Angular disparity - ",num2str(cfg_iRT_sessionID)," ",cfg_iRT_taskID," + ",seriesName];

  % subplot: rows, columns, index of selected plot
  #subplot (sub_p(1),sub_p(2),sub_p(3));
  [ax,h1,h2] = plotyy (x_duration , angDisp, x_duration, extra_series);
  title(plot_angDisp_title);
  xlabel("Task duration");
  ylabel (ax(1), "Angular disparity in degrees");
  ylabel (ax(2), seriesName);
  axis (ax(1), [0,Inf, 0,180] );
  axis (ax(2), 'autoy' ); % auto specifies the y-axis length
endfunction
% compute cummulative angular disparity in time
function series_cumulative = seriesSum(series)
  cumulativeSum = 0;
  series_cumulative = zeros(size(series));
  % 1st cumulative value is always zero.
  for i = 2:length(series_cumulative)
    cumulativeSum += abs (series(i-1) - series(i));
    series_cumulative(i) = cumulativeSum;
  endfor
endfunction
% compute cummulative smoothed angular disparity in time
function angDisp_cumulative_smooth = angDispSum2(angDisp, mw_radius)
  angDisp_smooth = angDisp;
  last_i = length(angDisp);
  if mw_radius < 1
    mw_radius = 1; %1 or higher for moving window radius
  endif
  for i = 1+mw_radius:last_i-mw_radius
    angDisp_smooth(i) = ( angDisp(i-1)+angDisp(i)+angDisp(i+1) ) /3;
  endfor
  cumulativeSum = 0;
  angDisp_cumulative_smooth = zeros(size(angDisp_smooth));
  % 1st cumulative value is always zero.
  for i = 2:length(angDisp_cumulative_smooth)
    cumulativeSum += abs (angDisp_smooth(i-1) - angDisp_smooth(i));
    angDisp_cumulative_smooth(i) = cumulativeSum;
  endfor
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
% draw atoms and gaze positions from both windows, between gaze_frames(a:b). frame is for the interactive model rotation state
function draw_atom_and_gaze (model_file_name, gaze_px, atom_px_tgt, atom_px_int, gaze_frames, frame = 1, cfg_iRT_sessionID, cfg_iRT_taskID)
  % frame is to select the interactive model orientation state
  %gaze_frames is the range of frame for the gaze squares to appear in the graph
  % read model file content
  [atom_count, atom_elem, atom_xyz] = get_xyz_data (strcat ("models/", model_file_name) );
  % build color and size vector for atoms
  atom_cor = generate_color_vector (atom_count, atom_xyz, atom_elem);
  figure();

  s1=subplot(1,2,1);
  %2D scatter of vertices/atoms
  scatter(atom_px_tgt(:,1),atom_px_tgt(:,2), atom_cor(:,1), atom_cor(:,2:4));
  axis ([300,900, 180,780] );
  set (gca (), "ydir", "reverse") %invert y axis so (0,0) is in top-left, like the screen pixels
  hold on;
  % scatter for gaze positions.
  scatter(gaze_px(gaze_frames,1),gaze_px(gaze_frames,2), 60, "g", "square");
  title ('Target window');
  xlabel("horizontal pixel position"); ylabel("vertical pixel position");
  hold off;

  s2= subplot(1,2,2);
  scatter(atom_px_int(:,1,frame),atom_px_int(:,2,frame), atom_cor(:,1), atom_cor(:,2:4));
  axis ([1140, 1740, 180,780] );
  set (s2, "ydir", "reverse") %invert y axis
  hold on;
  scatter(gaze_px(gaze_frames,1),gaze_px(gaze_frames,2), 60, "g", "square");
  set(s2, 'title', ['Interactive window at frame=',num2str(frame)]);
  xlabel("horizontal pixel position"); ylabel("vertical pixel position");
  legend("Vertices/Atoms", "Gaze");


  mainTitle = ["Atom & gaze screen projections,\n frames ", num2str(gaze_frames(1)),":",num2str(gaze_frames(length(gaze_frames)))," - ",num2str(cfg_iRT_sessionID)," ",cfg_iRT_taskID];
  S  = axes( 'visible', 'off', 'title', mainTitle );
  hold off;
  % sample call:
  %draw_atom_and_gaze (model_file_name, gaze_px, atom_px_tgt, atom_px_int, [1:30], time = 30, cfg_iRT_sessionID, cfg_iRT_taskID);
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
function atom_xyzRot_int = rotate_atom_xyz (task_data, atom_matrix_quat_cols, atom_xyz)
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
      atom_xyzRot_int(a,1:3,t) = (rot_vector(1:3,1:3,t)*atom_xyz(a,1:3)' )' ;
    endfor
  endfor
endfunction

% rotate atom_xyz atom matrix based on the target model quaternion from session_data
function atom_xyzRot_tgt = rotate_tgt_atom_xyz (session_data, session_row, cfg_atom_matrix_tgt_cols, atom_xyz)
  Q_tgt = cell2mat ( session_data(session_row,cfg_atom_matrix_tgt_cols) );
  for a = 1:size(atom_xyz,1); %all lines(atoms)
    atom_xyzRot_tgt(a,1:3) = ( rot_matrix (Q_tgt(1), Q_tgt(2), Q_tgt(3), Q_tgt(4)) * atom_xyz(a,1:3)' )' ;
  endfor
endfunction
% calculate, in pixels, the canvas center from the browser window
function canvas_center = get_cvs_center (session_data, row_index, cvs_pos_cols)
  canvas = cell2mat (session_data (row_index,cvs_pos_cols) ); %top, right, bottom, left side positions.
  canvas_center = [ (canvas(2) + canvas(4) )/2 , (canvas(1) + canvas(3) )/2 ];
endfunction
% compute pixel distance between gaze matrix and each atom of selected atom matrix
function [gaze_atomInt_dist, gaze_atomTgt_dist] = get_gaze_atom_dist (gaze_px, atom_px_int, atom_px_tgt)
%    printf("Calculating gaze-atom distance array.."); tic();
    %initialize array
    atom_count = size(atom_px_int,1);
    row_count = size(gaze_px,1);
    gaze_atom_dist = zeros ( row_count, atom_count );
    % compute matrices of screen position distances between gaze and each atom
    for a=1:atom_count %for each atom
      for t=1 : row_count %for each point in time
        %Formula: gaze_atom_dist = sqrt( (x-x')^2 + (y-y')^2 )
%        gaze_atom_dist(t,a) = sqrt ( (atom_px(a,1,t)-gaze_px(t,1))^2 + (atom_px(a,2,t)-gaze_px(t,2))^2 );
        gaze_atomInt_dist(t,a) = sqrt ( (atom_px_int(a,1,t)-gaze_px(t,1))^2 + (atom_px_int(a,2,t)-gaze_px(t,2))^2 );
        gaze_atomTgt_dist(t,a) = sqrt ( (atom_px_tgt(a,1)-gaze_px(t,1))^2 + (atom_px_tgt(a,2)-gaze_px(t,2))^2 );
      endfor
    endfor
%    printf(".array calculated.");
endfunction

% returns gaze_region (if gaze is inside canvas 1 or 2) from gaze pixel position and canvas extremities position in px [top right bottom left]
function gaze_region = fill_gaze_region (gaze_px, canvas_tgt, cfg_gaze_region_codeTgt, canvas_int, cfg_gaze_region_codeInt)
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
  %fill gaze_region if gaze is within canvas bounduaries
  for t=1 : frame_count
    % if gaze is inside target canvas..
    if ( (cvsLeft_r < gaze_px(t,1)) && (gaze_px(t,1) < cvsRight_r) && (cvsTop_r < gaze_px(t,2)) && (gaze_px(t,2) < cvsBottom_r) )
      gaze_region(t,1) = cfg_gaze_region_codeTgt;
    %if gaze is inside interactive canvas..
    elseif ( (cvsLeft_i < gaze_px(t,1)) && (gaze_px(t,1) < cvsRight_i) && (cvsTop_i < gaze_px(t,2)) && (gaze_px(t,2) < cvsBottom_i) )
      gaze_region(t,1) = cfg_gaze_region_codeInt;
    else
      gaze_region(t,1) = 0;
    endif
  endfor
endfunction
% stores nearest atom index for tgt and int. needs gaze_region and gaze_atomInt_dist (or Tgt)
function [gaze_nearest_atom_index_tgt,gaze_nearest_atom_index_int] = get_nearest_atom (gaze_region, gaze_atomTgt_dist, gaze_atomInt_dist)
  frame_count = size(gaze_region, 1);
  gaze_nearest_atom_index_tgt = zeros (frame_count,1);
  gaze_nearest_atom_index_int = zeros (frame_count,1);
  for t=1 : frame_count
    if gaze_region(t) == 1
      [gaze_nearest_atom_distance(t), gaze_nearest_atom_index_tgt(t)] =  min (gaze_atomTgt_dist(t,:));
    elseif gaze_region(t) == 2
      [gaze_nearest_atom_distance(t), gaze_nearest_atom_index_int(t)] =  min (gaze_atomInt_dist(t,:));
    endif
  endfor
endfunction
% compute the avg distance from all gaze_nearest_atoms inside the moving window to the center of them
function gaze_avg_dist = get_gaze_avg_dist (gaze_nearest_atom_index, atom_xyz, window_size) %DONT FORGET TO SET THE MODEL OF ATOMS (tgt or int)
  frame_count = size(gaze_nearest_atom_index, 1);
  center_xyz = zeros (frame_count,3);
  for t=1 : frame_count
    s = t-window_size; %the start of the moving window. Can't be lower than 1.
    if s < 1
      s=1;
    endif
    %get unique list of atom index [0,23,0,0,34,34]->[0,23,34]
    atom_list = unique ( gaze_nearest_atom_index(s:t) );
    % if the list has something other than [0]
    if ( length(atom_list) > 1 || atom_list(1)>0 )
      % if there is 0 with other indexes, remove 0 and keep computing [0,23,45,62]->[23,45,62]
      if (atom_list(1) == 0 )
        atom_list = atom_list(2:length(atom_list));
      endif
      center_xyz(t,1:3) = mean ( atom_xyz(atom_list,:) , 1) ; %compute center coordinate
      atom_list_dist = zeros (length(atom_list),1);
      % compute distance from center for each atom in list
      for i=1 : length (atom_list)
        atom_list_dist(i) = sqrt ( (atom_xyz(atom_list(i),1)-center_xyz(1))^2 + (atom_xyz(atom_list(i),2)-center_xyz(2))^2 + (atom_xyz(atom_list(i),3)-center_xyz(3))^2 );
      endfor
    elseif ( atom_list(1) == 0 && length(atom_list) == 1) %if there is only index 0, skip things. [0]
      center_xyz(t,1:3)  = [0,0,0];
      atom_list_dist = 0;
    endif
    gaze_avg_dist(t) = mean(atom_list_dist);
  endfor
endfunction
%plot the distance values of the gazemap matrix in scatter3 %debug
function plot_scatter3d_gazemap (gazemap_matrix_xxx)
  [gazemap_frames,gazemap_atoms] = size(gazemap_matrix_xxx);
  points_list = zeros (gazemap_atoms*gazemap_frames, 3);
  for i=1:gazemap_atoms
    for j=1:gazemap_frames
      testLine = gazemap_frames*(i-1) + j;
      points_list(testLine,1) = gazemap_matrix_xxx(j,i);
      points_list(testLine,2) = j;
      points_list(testLine,3) = i;
    endfor
  endfor
  figure ();
  scatter3 ( points_list(:,1), points_list(:,2), points_list(:,3) );
  frame_count = gazemap_frames;
  plot_title = ["3D Scatter of gazemap"];
    title (plot_title);
    xlabel("gazemap value"); ylabel("frame count"); zlabel("Atom index")
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

% returns jmol commands cell matrix for atoms transparency animation from gazemap_matrix_int (or _tgt). Needs to horzcat() for replay_int_jmol_script
function replay_transp_jmol_script = replay_transparency (frame_count, atom_count, heatmap_matrix)
  % building heatmap scale for animation
  for t=1:frame_count
    heatmap_mw_max(t) = max(heatmap_matrix(t,:));
    heatmap_mw_min(t) = min(heatmap_matrix(t,:));
    heatmap_mw_range(t) = heatmap_mw_max(t) - heatmap_mw_min(t);
    % compute scale from 1 (max) to nearly 0 for all atoms in each time frame
    for a=1:atom_count
      % 0.05 is minimal value to avoid division by 0
      heatscale_mw(t,a) = heatmap_matrix(t,a) / max ( heatmap_mw_max(t), 0.05 );
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
function writeOutput_gazemap (filename, script_matrix_tgt, script_matrix_int, script_line_tgt, script_line_int)
  tic(); printf(["Writing jmol commands in file '",filename,"' for gaze heatmap animation replay .."]);
  %open file pointer
  xls_gazemap_pointer = xlsopen (filename, true);

  [xls_gazemap_pointer] = oct2xls (script_matrix_int, xls_gazemap_pointer, "gaze replay int");
  [xls_gazemap_pointer] = oct2xls (script_matrix_tgt, xls_gazemap_pointer, "gaze replay tgt");
  [xls_gazemap_pointer] = oct2xls (script_line_int, xls_gazemap_pointer, "gaze frame int");
  [xls_gazemap_pointer] = oct2xls (script_line_tgt, xls_gazemap_pointer, "gaze frame tgt");

  %close file pointer
  [xls_gazemap_pointer] = xlsclose (xls_gazemap_pointer);
  printf(".OK! "); toc();
endfunction
% calculate gaze center from pondered gazemap. Compute pondered mean distance from previous calculated center. do for each frame in a window.
function [gazemap_center_xyz,gazemap_center_mean_disp] = compute_gazemap_center_dispersion (atom_xyz, gaze_region, gazemap_matrix_xxx)
  [frame_count, atom_count] = size(gazemap_matrix_xxx);
  gazemap_center_xyz = zeros(frame_count, 3);
  gazemap_center_mean_disp = zeros(frame_count, 1);
  frame_gaze_sum = zeros(frame_count, 1);
  for t = 1:frame_count
    % for each frame, obtain the mean contribution of all model atoms in the
    %matrix pondered by the gazemap_matrix (amount of gaze/attention received)
    % NOTE:gazemap_matrix_xxx is transposed to match the matrices dimensions in
    %the product (atoms as lines, xyz as columns).
    frame_gaze_sum(t) = sum ( gazemap_matrix_xxx( t,: ), 2 );
    if frame_gaze_sum(t) < 1 && t >= 2 && frame_gaze_sum(t-1) > frame_gaze_sum(t)
      %if the gaze center was "lost in the gaze" use previous value. Conditions:
      %1.all atom gaze values became extremely low in the moving window; 2.it is
      % not the 2 1st values; 3.the values are not decreasing.
      gazemap_center_xyz(t,1:3) =  gazemap_center_xyz(t-1,1:3);
    else
      gazemap_center_xyz(t,1:3) =  mean ( (gazemap_matrix_xxx( t,: )' .* atom_xyz(:,:) ) , 1) / sum ( gazemap_matrix_xxx( t,: )' );
    endif
    % compute distance from center for each atom in list
    gm(1,1:atom_count) = gazemap_matrix_xxx( t,: )';
    gazemap_center_dist = zeros(atom_count,1);
    for i=1 : atom_count
      gazemap_center_dist(i) = sqrt ( (gm(i).*atom_xyz(i,1)-gazemap_center_xyz(1))^2 + (gm(i).*atom_xyz(i,2)-gazemap_center_xyz(2))^2 + (gm(i).*atom_xyz(i,3)-gazemap_center_xyz(3))^2 );
    endfor
    gazemap_center_mean_disp(t) = mean(gazemap_center_dist);
  endfor
endfunction
% compute the translation speed of gazemap pondered center
function gazemap_center_speed = compute_gazemap_center_speed ( gazemap_center_xyz_xxx)
  time_frame = 0.1;
  gzm = gazemap_center_xyz_xxx;
  gazemap_center_speed = zeros (length ( gzm ), 1);
  for t = 2:length(gzm)
    gazemap_center_speed(t) = sqrt ( (gzm(t-1,1)-gzm(t,1))^2 + (gzm(t-1,2)-gzm(t,2))^2 + (gzm(t-1,3)-gzm(t,3))^2 )/time_frame;
  endfor
endfunction
%compute distance between gaze centers
function gazemap_centers_delta = get_gaze_center_deltas (gazemap_center_xyz_tgt, gazemap_center_xyz_int)
  frames = length(gazemap_center_xyz_tgt);
  int = gazemap_center_xyz_int;
  tgt = gazemap_center_xyz_tgt;
  gazemap_centers_delta = zeros(frames, 1);
  for i=1:frames
    gazemap_centers_delta(i) = sqrt ( (tgt(i,1)-int(i,1))^2 + (tgt(i,2)-int(i,2))^2 + (tgt(i,3)-int(i,3))^2  );
  endfor
endfunction
% scatter in 2x2 grid for data split in stage (seek/adjustment) and canvas (target/interactive)
function scatterplot_2x2 (swap_frame, x_tgt,x_int, y_tgt,y_int, x_label,y_label, color_data = 0)
  frame_count = length(x_tgt);
  mrk_size = 8; %marker size, default is 36
  if color_data == 0
    mrk_color = [0:0.1:length(x_tgt)]'; %marker color scheme
  else
    mrk_color = color_data;
  endif
  %colorbar adds the "color legend"
  figure;
  subplot(2,2,1);
    h1 = scatter ( x_tgt(1:swap_frame) , y_tgt(1:swap_frame) , mrk_size, mrk_color(1:swap_frame));
    title ("Target, seek stage");
    xlabel(x_label); ylabel(y_label);
    colorbar("EastOutside");
  subplot(2,2,2);
    h1 = scatter ( x_tgt(swap_frame:frame_count) , y_tgt(swap_frame:frame_count) , mrk_size, mrk_color(swap_frame:frame_count));
    title ("Target, adjustments stage");
    xlabel(x_label); ylabel(y_label);
    colorbar("EastOutside");
  subplot(2,2,3);
    h3 = scatter ( x_int(1:swap_frame) , y_int(1:swap_frame) , mrk_size, mrk_color(1:swap_frame) );
    title ("Interactive, seek stage");
    xlabel(x_label); ylabel(y_label);
    colorbar("EastOutside");
  subplot(2,2,4);
    h4 = scatter ( x_int(swap_frame:frame_count) , y_int(swap_frame:frame_count) , mrk_size, mrk_color(swap_frame:frame_count));
    title ("Interactive, adjustments stage");
    xlabel(x_label); ylabel(y_label);
    colorbar("EastOutside");

endfunction
%scatter in 2x3 grid for data split in stage (seek/adjustment/complete) and canvas (target/interactive)
function scatterplot_2x3 (swap_frame, x_tgt,x_int, y_tgt,y_int, x_label,y_label, color_data = 0)
  frame_count = length(x_tgt);
  mrk_size = 8; %marker size, default is 36
  if color_data == 0
    mrk_color = [0:0.1:length(x_tgt)]'; %marker color scheme
  else
    mrk_color = color_data;
  endif
  %colorbar adds the "color legend"
  figure;
  subplot(2,3,1);
    h1 = scatter ( x_tgt(1:swap_frame) , y_tgt(1:swap_frame) , mrk_size, mrk_color(1:swap_frame));
    title ("Target, seek stage");
    xlabel(x_label); ylabel(y_label);
    colorbar("EastOutside");
  subplot(2,3,2);
    h2 = scatter ( x_tgt(swap_frame:frame_count) , y_tgt(swap_frame:frame_count) , mrk_size, mrk_color(swap_frame:frame_count));
    title ("Target, adjustments stage");
    xlabel(x_label); ylabel(y_label);
    colorbar("EastOutside");
  subplot(2,3,3);
    h3 = scatter ( x_tgt(1:frame_count) , y_tgt(1:frame_count) , mrk_size, mrk_color(1:frame_count));
    title ("Target, complete");
    xlabel(x_label); ylabel(y_label);
    colorbar("EastOutside");
  subplot(2,3,4);
    h4 = scatter ( x_int(1:swap_frame) , y_int(1:swap_frame) , mrk_size, mrk_color(1:swap_frame) );
    title ("Interactive, seek stage");
    xlabel(x_label); ylabel(y_label);
    colorbar("EastOutside");
  subplot(2,3,5);
    h5 = scatter ( x_int(swap_frame:frame_count) , y_int(swap_frame:frame_count) , mrk_size, mrk_color(swap_frame:frame_count));
    title ("Interactive, adjustments stage");
    xlabel(x_label); ylabel(y_label);
    colorbar("EastOutside");
  subplot(2,3,6);
    h6 = scatter ( x_int(1:frame_count) , y_int(1:frame_count) , mrk_size, mrk_color(1:frame_count));
    title ("Interactive, complete");
    xlabel(x_label); ylabel(y_label);
    colorbar("EastOutside");
endfunction
%scatter in 1x3 grid
function scatterplot_1x3 (swap_frame, x_var, y_var, x_label,y_label, color_data = 0)
  frame_count = length(x_var);
  mrk_size = 8; %marker size, default is 36
  if color_data == 0
    mrk_color = [0:0.1:length(x_var)]'; %marker color scheme
  else
    mrk_color = color_data;
  endif
  %colorbar adds the "color legend"
  figure;
  subplot(1,3,1);
    h1 = scatter ( x_var(1:swap_frame) , y_var(1:swap_frame) , mrk_size, mrk_color(1:swap_frame));
    title ("seek stage");
    xlabel(x_label); ylabel(y_label);
    colorbar("EastOutside");
  subplot(1,3,2);
    h2 = scatter ( x_var(swap_frame:frame_count) , y_var(swap_frame:frame_count) , mrk_size, mrk_color(swap_frame:frame_count));
    title ("adjustments stage");
    xlabel(x_label); ylabel(y_label);
    colorbar("EastOutside");
  subplot(1,3,3);
    h3 = scatter ( x_var(1:frame_count) , y_var(1:frame_count) , mrk_size, mrk_color(1:frame_count));
    title ("complete");
    xlabel(x_label); ylabel(y_label);
    colorbar("EastOutside");

endfunction
%scatter in 2x3 grid for data split in stage (seek/adjustment/complete) but same x scale
function scatterplot_2x3_2ys (swap_frame, x_tgt,x_int, y_tgt,y_int, x_label,y_label, y_labelB, color_data = 0)
  frame_count = length(x_tgt);
  mrk_size = 8; %marker size, default is 36
  if color_data == 0
    mrk_color = [0:0.1:length(x_tgt)]'; %marker color scheme
  else
    mrk_color = color_data;
  endif
  %colorbar adds the "color legend"
  figure;
  subplot(2,3,1);
    h1 = scatter ( x_tgt(1:swap_frame) , y_tgt(1:swap_frame) , mrk_size, mrk_color(1:swap_frame));
    title ("seek stage");
    xlabel(x_label); ylabel(y_label);
    colorbar("EastOutside");
  subplot(2,3,2);
    h2 = scatter ( x_tgt(swap_frame:frame_count) , y_tgt(swap_frame:frame_count) , mrk_size, mrk_color(swap_frame:frame_count));
    title ("adjustments stage");
    xlabel(x_label); ylabel(y_label);
    colorbar("EastOutside");
  subplot(2,3,3);
    h3 = scatter ( x_tgt(1:frame_count) , y_tgt(1:frame_count) , mrk_size, mrk_color(1:frame_count));
    title ("complete");
    xlabel(x_label); ylabel(y_label);
    colorbar("EastOutside");
  subplot(2,3,4);
    h4 = scatter ( x_int(1:swap_frame) , y_int(1:swap_frame) , mrk_size, mrk_color(1:swap_frame) );
    title ("seek stage");
    xlabel(x_label); ylabel(y_labelB);
    colorbar("EastOutside");
  subplot(2,3,5);
    h5 = scatter ( x_int(swap_frame:frame_count) , y_int(swap_frame:frame_count) , mrk_size, mrk_color(swap_frame:frame_count));
    title ("adjustments stage");
    xlabel(x_label); ylabel(y_labelB);
    colorbar("EastOutside");
  subplot(2,3,6);
    h6 = scatter ( x_int(1:frame_count) , y_int(1:frame_count) , mrk_size, mrk_color(1:frame_count));
    title ("complete");
    xlabel(x_label); ylabel(y_labelB);
    colorbar("EastOutside");
endfunction
% discretize data. bin_size is the numeric size of each discrete "step" (ex: 0,15,30,45 for a size of 15).
function data_discrete = discretize_data (data, bin_size)
  data_discrete = data /bin_size;
  data_discrete = fix(data_discrete );
  data_discrete = data_discrete *bin_size;
endfunction

% moving average of gaze_region plot. mw_size is amount of points including center
function gaze_region_mw = get_gaze_region_mw (mw_size, gaze_region, both_sides=0)
  frame_count = length(gaze_region);
  gaze_region_mw = zeros(frame_count,1);
  gs_arr = gaze_region;
  tgt_id = 1;
  int_id = 2;
  % rewrite gaze_region in gs_arr as -1 to 1, and empty when out of the canvases (because empty is skipped in mean func.)
  for i=1:frame_count
    if gaze_region(i) == tgt_id
      gs_arr(i) = -1;
    elseif gaze_region(i) == int_id
      gs_arr(i) = 1;
    else
      gs_arr(i) = 0;
    endif
  endfor
%  figure;plot(gs_arr); %DEBUG. should vary between -1 and 1
  % set moving average loop
  for i=mw_size:frame_count
    if i < mw_size
      gaze_region_mw(i) = mean ( gs_arr( 1:i ) );
    else
      gaze_region_mw(i) = mean ( gs_arr( 1+i-mw_size:i ) );
    endif
  endfor
endfunction
%cumulative and percent plot of gaze_region
function [gaze_region_sum,gaze_region_percent ]= get_gaze_region_sum (gaze_region, cfg_t,cfg_iRT_sessionID, cfg_iRT_taskID)
  frame_count = length(gaze_region);
  gaze_region_sum = zeros ( frame_count, 3);
  gaze_region_percent = zeros ( frame_count, 3);
  %setup the initial value
  switch ( gaze_region(1) )
    case {0} %outside
      gaze_region_sum(1,3) = cfg_t;
    case {1} %target
      gaze_region_sum(1,1) = cfg_t;
    case {2} %interactive
      gaze_region_sum(1,2) = cfg_t;
  endswitch
  %sum all remaining values
  for i=2:frame_count
    switch ( gaze_region(i) )
      case {0} %outside
        gaze_region_sum(i,3) = gaze_region_sum(i-1,3) + cfg_t;
        gaze_region_sum(i,1) = gaze_region_sum(i-1,1);
        gaze_region_sum(i,2) = gaze_region_sum(i-1,2);
      case {1} %target
        gaze_region_sum(i,3) = gaze_region_sum(i-1,3);
        gaze_region_sum(i,1) = gaze_region_sum(i-1,1) + cfg_t;
        gaze_region_sum(i,2) = gaze_region_sum(i-1,2);
      case {2} %interactive
        gaze_region_sum(i,3) = gaze_region_sum(i-1,3);
        gaze_region_sum(i,1) = gaze_region_sum(i-1,1);
        gaze_region_sum(i,2) = gaze_region_sum(i-1,2) + cfg_t;
    endswitch
  endfor
  plot_temporal_data (gaze_region_sum, "Cumulative regional gaze", cfg_iRT_sessionID, cfg_iRT_taskID);
  legend("Target model", "Interactive Model", "Outside");

  %{
  %we did not found the following ratio usefull

  %percentage value (ratio) of the gaze sum
  for i=1:frame_count
    gaze_region_percent(i,:) = gaze_region_sum(i,:)./(i*cfg_t);
  endfor
  plot_temporal_data (gaze_region_percent, "Percent regional gaze", cfg_iRT_sessionID, cfg_iRT_taskID)
  legend("Target model", "Interactive Model", "Outside");

  %percentage value (ratio) of the gaze sum without including the "outside" gaze region
  gaze_region_sum_clean = sum(gaze_region_sum(:,1:2), 2);
  for i=1:frame_count
    gaze_region_percent_clean(i,1:2) = gaze_region_sum(i,1:2)./(gaze_region_sum_clean(i));
    if gaze_region_sum_clean(i) == 0
      gaze_region_percent_clean(i,1:2)=0.5; %or 0.5
    endif
  endfor
  plot_temporal_data (gaze_region_percent_clean, "Percent regional gaze", cfg_iRT_sessionID, cfg_iRT_taskID);
  legend("Target model", "Interactive Model");
  %}

  %plot of gaze_region as ref x int
  plot_data_OLD (gaze_region_sum(:,2),"Gaze Int", gaze_region_sum(:,1),"Gaze Ref", cfg_iRT_sessionID, cfg_iRT_taskID,0,0);

endfunction


% compute gaze deslocation (how much pixels it moves per time-frame
function [gaze_px_movement,gaze_px_dx,gaze_px_dy] = compute_gaze_px_movement (gaze_px)
  frame_count = length(gaze_px);
  gaze_px_dx = gaze_px_dy = gaze_px_movement = zeros(frame_count,1);
  gaze_px_dx(2:frame_count,1) = gaze_px(2:frame_count, 1) - gaze_px(1:frame_count-1, 1);
  gaze_px_dy(2:frame_count,1) = gaze_px(2:frame_count, 2) - gaze_px(1:frame_count-1, 2);
  gaze_px_movement(2:frame_count) = norm ( [ gaze_px_dx(2:frame_count), gaze_px_dy(2:frame_count) ], 2, "rows");
  for i=2:frame_count
    gaze_px_movement(i) = sqrt ( gaze_px_dx(i)^2 + gaze_px_dy(i)^2 );
  endfor
endfunction
%TODO: write a function for plotting available graphs
%==========================



%SCRIPTS
%--------------
% 1.Temporal Data Input
printf("\n STEP 1: Temporal Data Input\n");
% data checks for the slowest functions! If the variable does not exist and is
%configured to not get new data, it changes config to get new data
if and ( exist('raw_iRT_data', 'var') == 0 , cfg_iRT_input == false)
  warning("The iRT data source is missing! Changed 'cfg_iRT_input' to true for this execution.\n");
  cfg_iRT_input = true;
endif
if and ( exist('raw_eyeT_data', 'var') == 0 , cfg_eyeT_input == false)
  warning("The eyeTracking data source is missing! Changed 'cfg_eyeT_input' to true for this execution.\n");
  cfg_eyeT_input = true;
endif

% INPUT ROTATIONAL DATA
% iRT_data and session_data input (always computes session_row, Q_tgt, Q and frame_count)
if cfg_iRT_input == true
  [session_data,raw_iRT_data] = iRT2oct (cfg_iRT_input_filename);
else
  disp("Skipping iRT file read.");
endif

% filter values from specified sessionID (cfg_iRT_sessionID) and taskID (cfg_iRT_taskID) inside session_data
if cfg_iRT_filter == true
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
    plot_temporal_data (angDisp, "Angular Disparity in degrees", cfg_iRT_sessionID, cfg_iRT_taskID, [0,180])
    legend('Angular Disparity');
  endif
endif

% INPUT GAZE DATA
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
% ----------------------
% 2.Data Synchronization
printf("\n STEP 2: Data Synchronization\n");
% synchronization/merge of iRT and eyeTracking data. Headers included. 1st column values of matrices is used for synchronization
if cfg_data_merge == true
  % merge eyeT data into iRT data
  task_data = merge_data(iRT_data, eyeT_data, cfg_eyeT_cols);
  % merge headers
  merged_header_data = horzcat(iRT_header_data, "angDisp", raw_eyeT_header_data(cfg_eyeT_cols));

  % WRITE output file
  if cfg_write_merge_output == true
    writeOutput_merge_fname = ["output merge ", num2str(cfg_iRT_sessionID), " ", cfg_iRT_taskID, ".xlsx"];
    writeOutput_merged (["output merge ", num2str(cfg_iRT_sessionID), " ", cfg_iRT_taskID, ".xlsx"], merged_header_data, task_data, session_data, session_row);
  endif
endif
%------------------
% 3.Data Processing
%TODO: find basic data variables and compute them first. example: gaze_px
printf("\n STEP 3: Data Processing\n");
% Build interaction replay animation from temporal quaternion (t x 4 array)
if cfg_replay_animation == true
  cfg_replay_animation_filename = ["output jmol console ",num2str(cfg_iRT_sessionID)," ",cfg_iRT_taskID,".xlsx"];
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
  tic(); printf("Writing file with jmol commands for replay animation of rotations..");
  %open file pointer
  xls_replay_pointer = xlsopen (cfg_replay_animation_filename, true);

  [xls_replay_pointer] = oct2xls (replay_int_jmol_script, xls_replay_pointer, 'rotaion replay');

  %close file pointer
  [xls_replay_pointer] = xlsclose (xls_replay_pointer);
  printf(".OK!"); toc();
endif


% CUMULATIVE series of Angular Disparity.
if cfg_angDisp_cumulative == true
  angDisp_cumulative = seriesSum(angDisp);
  if cfg_plot_angDisp_cumulative == true
    plot_angDisp_yy_other (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, angDisp_cumulative, 'Cumulative angDisp');
  endif
  % discrete cumulative series of angular disparity
  angDisp_discrete = discretize_data (angDisp, 15);
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
    title ("3D colored Scatter of vertices");
    axis ("equal");
    xlabel("x"); ylabel("y"); zlabel("z");
  endif
endif

% compute rotated atom matrices (temporal and target) from xyz coordinates rotated acording to quaternion data in chosen session
if cfg_atom_matrix == true
  % create temporal atom matrix of the interactive model
  atom_xyzRot_int = rotate_atom_xyz (task_data, cfg_atom_matrix_quat_cols, atom_xyz);
  % create atom(vertices) matrix of the target model
  atom_xyzRot_tgt = rotate_tgt_atom_xyz (session_data, session_row, cfg_atom_matrix_tgt_cols, atom_xyz);

  % compute xy pixel coordinate of canvas center (both target and interactive)
  cvsInt_center = get_cvs_center (session_data, session_row, cfg_gaze_cvsInt_cols);
  cvsTgt_center = get_cvs_center (session_data, session_row, cfg_gaze_cvsTgt_cols);
  % get pixel to angstrom ratio of chosen row in session_data
  pxAngs_rate = session_data{session_row, cfg_gaze_pxAngs_rate_col};
  if ischar(pxAngs_rate)
    pxAngs_rate = str2num(pxAngs_rate); % string -> number
  endif
  % get pixel ratio on screen of chosen row in session_data
  pxRatio = session_data{session_row, cfg_gaze_pxRatio_col};
  if ischar(pxRatio)
    pxRatio = str2num(pxRatio);
  endif

  % get px coordinates of atoms xy projection (temporal/interactive and static/target).
  %NOTE: (0,0) coodinate is the top-left screen corner, but in Jmol it is the bottom-left corner.
  % Hence, Jmol y coordinates must be inverted before summing with the canvas center.
  % pxAngs_rate is the ration between Angstrons (Jmol measure) and pixels on screen,
  % and pxRatio is the ratio of pixels on screen (from the windows resource used mostly in notebooks).
  atom_px_int = (atom_xyzRot_int(:,1:2,:) .* (pxAngs_rate * pxRatio) );
  atom_px_int(:,2,:) = -atom_px_int(:,2,:);
  atom_px_int = atom_px_int + cvsInt_center;
  atom_px_tgt = (atom_xyzRot_tgt(:,1:2) .* (pxAngs_rate * pxRatio) );
  atom_px_tgt(:,2) = -atom_px_tgt(:,2);
  atom_px_tgt = atom_px_tgt + cvsTgt_center;
  %plot interactive model atoms/vertices in 2D at given time set as cfg_atom_matrix_plot_t
  if and ( cfg_atom_matrix_plot == true, cfg_atom_matrix_plot_t > size(atom_xyzRot_int,3) )
    warning("The chosen frame index %i is outside the matrix range. Choose a reasonable frame index");
  elseif cfg_atom_matrix_plot == true
    atom_cor = generate_color_vector (atom_count, atom_xyz, atom_elem);
    figure (3);
    scatter (atom_xyzRot_int(:,1,cfg_atom_matrix_plot_t), atom_xyzRot_int(:,2,cfg_atom_matrix_plot_t), atom_cor(:,1), atom_cor(:,2:4));
    title ( strcat ("2D Scatter of vertices from interactive object at time=", num2str(cfg_atom_matrix_plot_t) ) );
    axis ("equal");
    xlabel("x"); ylabel("y");
  endif
  %plot target model atoms/vertices in 2D
  if cfg_atom_matrix_tgt_plot == true
    atom_cor = generate_color_vector (atom_count, atom_xyz, atom_elem);
    figure (4);
    scatter (atom_xyzRot_tgt(:,1), atom_xyzRot_tgt(:,2), atom_cor(:,1), atom_cor(:,2:4));
    title ("2D Scatter of vertices in target model");
    axis ("equal");
    xlabel("x"); ylabel("y");
  endif

endif
% plot angular disparity with PUPIL DIAMETER data (TODO: include pupil computation inside other set of gaze computations?)
if cfg_plot_angDisp_yy_pupil == true
  gaze_pupil_diameter_mean = task_data(:,13)+task_data(:,14)/2;
  plot_angDisp_yy_other (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, gaze_pupil_diameter_mean, "pupil data");
endif
% compute matrices of distances on screen between gaze and each atom
if cfg_gaze_dist_matrix == true
  gaze_px = task_data (:,cfg_gaze_cols);
  [gaze_atomInt_dist, gaze_atomTgt_dist] = get_gaze_atom_dist (gaze_px, atom_px_int, atom_px_tgt);
endif
% compute 'gaze_region' in relation to canva :  1= tgt; 2= int; 0= outside
if cfg_gaze_region_array == true
  canvas_int = cell2mat (session_data (session_row,cfg_gaze_cvsInt_cols) ); % top, right, bottom, left side positions.
  canvas_tgt = cell2mat (session_data (session_row,cfg_gaze_cvsTgt_cols) ); % top, right, bottom, left side positions.
  gaze_region = zeros ( size(gaze_px,1), 1); % n x 1
  gaze_region(:,1) = fill_gaze_region (gaze_px, canvas_tgt, cfg_gaze_region_codeTgt, canvas_int, cfg_gaze_region_codeInt);
  gaze_region_mw = get_gaze_region_mw (11, gaze_region);

  % plot angular disparity colored with gaze_region
  if cfg_plot_angDisp_gaze_region == true
    plot_angDisp_colored_bg (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, gaze_region);
  endif
  if cfg_plot_gaze_region_mw == true
    plot_angDisp_yy_other (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, gaze_region_mw, "Gaze region (moving-window)");
  endif
endif
% Calculate 3D gaze mapping using bigaussian for distance between gaze position and each atom on screen at each time, with normal (another gaussian) decay in time
% TODO: CHECK FUNCTIONS AND VARIABLES (GAZE_DIST ETC)
if cfg_gaze_map == true
  printf("Calculating transparency gradient for replay animation (may take a while):\n");tic();
  %declaring variables
  gazemap_matrix_int = gazemap_matrix_tgt = zeros (atom_count,1); #(atoms,1)
  gaze_quad_dist_tgt = gaze_quad_dist_int = zeros (frame_count,2,atom_count);  #(time,1:2,atoms)
  % porque (t,1:2,a) ? não deveria ser (t,a)?
%  gaze_quad_dist_tgt = gaze_quad_dist_int = zeros (frame_count,atom_count);  #(time,1:2,atoms)
  exp_gaze_dist_tgt = exp_gaze_dist_int = zeros (frame_count, atom_count); #(time,atoms)
  %Gaussian formula remainder:
  %integral of ( exp( - ( (x(t)-cx)^2 + (y(t)-cy)^2)/ (2*cfg_gauss_wdt_screen^2)) dt)
  %integral of ( exp( - ( (ax(t)-gx)^2 + (ay(t)-gy)^2)/ (2*cfg_gauss_wdt_screen^2)) dt)

  % Compute matrix of euclidian distance (gaze_quad_dist_tgt and gaze_quad_dist_int)
  %between gaze and atom, for all atoms in all time frames. ( (ax-gx)^2+(ay-gy)^2 )
  for a=1 : atom_count
    for t=1 : frame_count
      gaze_quad_dist_tgt(t,a) = [ ( atom_px_tgt(a,1)   - gaze_px(t,1) ).^2 + ( atom_px_tgt(a,2)   - gaze_px(t,2) ).^2 ];
      gaze_quad_dist_int(t,a) = [ ( atom_px_int(a,1,t) - gaze_px(t,1) ).^2 + ( atom_px_int(a,2,t) - gaze_px(t,2) ).^2 ];
    endfor
  endfor
  % fill bi-gaussian values for each atom (gaze distance gaussian)
  %(exp (- (gaze_dist)/(2*cfg_gauss_wdt_screen^2)
  for a=1 : atom_count
    exp_gaze_dist_tgt(:,a) = exp ( - gaze_quad_dist_tgt(:,a) / (2*cfg_gauss_wdt_screen^2) ) ;
    exp_gaze_dist_int(:,a) = exp ( - gaze_quad_dist_int(:,a) / (2*cfg_gauss_wdt_screen^2) ) ;
  endfor

  #defining visual short-term memory gaussian
  clear vstm_gaussian;
  for t=1 : 4*cfg_gauss_wdt_time
    vstm_gaussian(t) = exp ( - ( (4*cfg_gauss_wdt_time-t)^2 / (2*(cfg_gauss_wdt_time^2)) ) );
  endfor
  % plot(vstm_gaussian) % debug
  vstm_len = length(vstm_gaussian);

  % sum of gaze density/gaze mapping. For each atom, sums all values inside the bigaussian ("integrate"),
  % apply the memory gaussian (past time/memory degradation gaussian) and register
  % TODO : MAKE A FUNCTION FOR THIS!
  for t=1 : frame_count
    gazemap_sum_tgt = 0;
    gazemap_sum_int = 0;
    ti = 1 + max([t-4*cfg_gauss_wdt_time, 0]);
    for frame=ti : t
      gazemap_sum_tgt += ( exp_gaze_dist_tgt(frame,1:atom_count) .* vstm_gaussian(vstm_len-(t-frame)) );
      gazemap_sum_int += ( exp_gaze_dist_int(frame,1:atom_count) .* vstm_gaussian(vstm_len-(t-frame)) );
    endfor
    gazemap_matrix_tgt(t,1:atom_count) = gazemap_sum_tgt;
    gazemap_matrix_int(t,1:atom_count) = gazemap_sum_int;
  endfor

  % cumulative gaze mapping of the entire process, with no decay in time.
  gazemap_single_tgt(1:atom_count) = sum ( exp_gaze_dist_tgt(:,1:atom_count) );
  gazemap_single_int(1:atom_count) = sum ( exp_gaze_dist_int(:,1:atom_count) );

  % Build map of rotations with gaze transparency, if both gazemap and raplay animation are configured as true
  if cfg_replay_animation == true
    % building the list of Jmol console commands:

    % ..for target Object replay
    replay_tgt_jmol_script = repmat ({"delay 0.1;"}, frame_count, 1); %declaring with "delay 0.1"
    Q_tgt = cell2mat ( session_data (session_row,cfg_atom_matrix_tgt_cols) ); %"initial" Quaternion coordinate for target model
    replay_tgt_jmol_script{1} = ["moveto 0 QUATERNION {", num2str(Q_tgt(1:4)),"};"]; % 1st value as a moveto
    jmol_script_gazemap_replay_tgt = horzcat ( replay_tgt_jmol_script, replay_transparency (frame_count, atom_count, gazemap_matrix_tgt) );

    % ..for interactive Object replay
    %replay_int_jmol_script was already calculated in the script
    jmol_script_gazemap_replay_int = horzcat ( replay_int_jmol_script, replay_transparency (frame_count, atom_count, gazemap_matrix_int) );

    % ..for a single frame with each object (target and interactive)
    jmol_script_gazemap_single_tgt = horzcat ( replay_tgt_jmol_script(1), replay_transparency (1, atom_count, gazemap_single_tgt) );
    jmol_script_gazemap_single_int = horzcat ( strcat ("moveto 0.0 QUATERNION {", num2str ( Q(frame_count) ) , "};"), replay_transparency (1, atom_count, gazemap_single_int) );

    writeOutput_gazemap (cfg_replay_animation_filename, jmol_script_gazemap_replay_tgt, jmol_script_gazemap_replay_int, jmol_script_gazemap_single_tgt, jmol_script_gazemap_single_int);
  endif
endif
% compute the avg. dist. in time between center and recently gazed atoms (NO gazemap matrix)
if cfg_gaze_avg_dist == true
  %get nearest atom from gaze of the gazed canva.
  [gaze_nearest_atom_index_tgt,gaze_nearest_atom_index_int] = get_nearest_atom (gaze_region, gaze_atomTgt_dist, gaze_atomInt_dist);
  %Then, compute the avg dist in time of the space being formed between these gazed atoms.
  gaze_avg_dist_tgt = get_gaze_avg_dist (gaze_nearest_atom_index_tgt, atom_xyz, cfg_gaze_window_size);
  gaze_avg_dist_int = get_gaze_avg_dist (gaze_nearest_atom_index_int, atom_xyz, cfg_gaze_window_size);
  figure;
  plot_angDisp_yy_other (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, gaze_avg_dist_tgt, 'mean distance of tgt atoms', [3,1]);
  plot_angDisp_colored_bg (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, gaze_region, [3,2]);
  plot_angDisp_yy_other (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, gaze_avg_dist_int, 'mean distance of int atoms', [3,3]);
endif
%compute gaze dispersion (mean distance of each recently gazed atom from calculated center)
if cfg_gaze_dispersion == true
  [gazemap_center_xyz_tgt, gazemap_center_dispersion_tgt] = compute_gazemap_center_dispersion (atom_xyz, gaze_region, gazemap_matrix_tgt);
  [gazemap_center_xyz_int, gazemap_center_dispersion_int] = compute_gazemap_center_dispersion (atom_xyz, gaze_region, gazemap_matrix_int);
  if cfg_plot_gaze_dispersion == true;
    figure;
    plot_angDisp_yy_other (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, gazemap_center_dispersion_tgt, 'TGT Gaze Dispersion',[3,1]);
    plot_angDisp_colored_bg (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, gaze_region, [3,2]);
    plot_angDisp_yy_other (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, gazemap_center_dispersion_int, 'INT Gaze Dispersion',[3,3]);
  endif
  %compute distance between gaze centers
  if cfg_gaze_centers_delta == true
    gazemap_centers_delta = get_gaze_center_deltas (gazemap_center_xyz_tgt, gazemap_center_xyz_int);
    if cfg_plot_gaze_centers_delta == true
      plot_angDisp_yy_other (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, gazemap_centers_delta, 'Gaze centers distance(tgt-int)');
    endif
  endif
  %Compute gaze center speed (deslocation / time)
  if cfg_gaze_dispersion_center_speed == true
    gazemap_center_speed_tgt = compute_gazemap_center_speed ( gazemap_center_xyz_tgt);
    gazemap_center_speed_int = compute_gazemap_center_speed ( gazemap_center_xyz_int);
    if cfg_plot_gaze_dispersion_center_speed == true
      figure;
      plot_angDisp_yy_other (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, gazemap_center_speed_tgt, 'TGT Gaze center speed', [3,1]);
      plot_angDisp_yy_other (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, gazemap_center_speed_int, 'INT Gaze center speed', [3,3]);
      plot_angDisp_yy_other (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, [seriesSum(gazemap_center_speed_tgt),seriesSum(gazemap_center_speed_int)], ['center spd TGT x INT'], [3,2]);
      legend_labels = {'Angular Disparity', 'Target model cumulative shift', 'Interactive Model cumulative shift'}; % is this working?
      legend ('Angular Disparity', 'Target model cumulative shift', 'Interactive model cumulative shift');
    endif
    % TODO : make second y axis range be the same, or change color of each line.
  endif
endif

% find the frame where exploration ends and adjustments begins. (TODO)
if cfg_finding_stage_change == true
  % TODO : remove functions inside condition when functions are giving expected results
  % compute simple 2-stages division in angular disparity plot (exploration X adjustments). It uses optim pkg
  function R_squared_val = compute_r_squared (series)
    frame_count = length(series);
    x = 0.1*(1:frame_count);
    y = series';

    % Fit a line (polynomial of degree 1) to the data
    [p, s] = polyfit(x, y, 1);

    % Calculate predicted y values using the fitted line
    y_pred = polyval(p, x);

    % Calculate total sum of squares (SST)
    y_mean = mean(y);
    SST = sum((y - y_mean).^2);

    % Calculate sum of squares of residuals (SSE)
    SSE = sum((y - y_pred).^2);

    % Calculate R-squared
    R_squared_val = 1 - (SSE / SST);
  endfunction

  % find the frame that best separates the angular disparity in two lines or stages (exploration X adjustments).
  function [R_sq_frame,R_sq_array] = find_aha_moment (angDisp_cumulative)
    %Returns int and array of the product of R_squared_val values through each chosen frame
    if length(angDisp_cumulative) < 7
      error("The array used is smaller than 7 points. This may return unexpected results. Are you sure this is the right data?\n");
    endif
    frame_count = size(angDisp_cumulative,1);
    safe = 4; %safety range, to avoid calculating the R-squared of arrays too small
    R_sq_frame = [0,0];
    R_sq_array = zeros(frame_count, 3); % for debug plots
    current_R_squared_val = 0;
    for frame = 1+safe : frame_count-safe
      R_stage1 = compute_r_squared (angDisp_cumulative(1:frame));
      R_stage2 = compute_r_squared (angDisp_cumulative(frame:frame_count));
      if R_stage1 > 0 && R_stage2 > 0
        R_sq_array(frame,1) = R_stage1;
        R_sq_array(frame,2) = R_stage2;
        R_sq_array(frame,3) = sqrt(R_stage1*R_stage2);
        if R_stage1*R_stage2 > current_R_squared_val
          current_R_squared_val = sqrt(R_stage1*R_stage2);
          R_sq_frame = [frame,current_R_squared_val];
        endif
      else
        R_sq_array(frame,1:3) = 0;
      endif
    endfor
  endfunction


  [R_sq_frame,R_sq_array] = find_aha_moment (angDisp_cumulative);
  figure; subplot(4,1,1); plot(R_sq_array(:,1));subplot(4,1,2); plot(R_sq_array(:,2));subplot(4,1,3); plot(R_sq_array(:,3)); plot_angDisp_colored_bg (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, gaze_region, sub= [4,4]);
endif

% visualize vertices and gaze xy positions
if cfg_draw_atom_and_gaze == true
  draw_atom_and_gaze (model_file_name, gaze_px, atom_px_tgt, atom_px_int, [1:30], time = 30, cfg_iRT_sessionID, cfg_iRT_taskID);
endif
% num. of gazed atoms in time (atom/gaze)
cfg_gaze_radius = 50
gazed_atoms = gazed_atoms_Int = gazed_atoms_Tgt = zeros(frame_count,1);
for t=1:frame_count
  gazed_atoms_Tgt (t,1) = sum ( gaze_atomTgt_dist(t,:) <= cfg_gaze_radius);
  gazed_atoms_Int (t,1) = sum ( gaze_atomInt_dist(t,:) <= cfg_gaze_radius);
  gazed_atoms (t,1) = gazed_atoms_Int (t,1) + gazed_atoms_Tgt (t,1);
endfor
  figure;
  plot_temporal_data (gazed_atoms_Tgt, "Target atoms", cfg_iRT_sessionID, cfg_iRT_taskID, axis_sizes=0, ytick_val=0, fig_num=[3,1,1]);
  plot_temporal_data (gazed_atoms_Int, "Interactive atoms", cfg_iRT_sessionID, cfg_iRT_taskID, axis_sizes=0, ytick_val=0, fig_num=[3,1,2]);
  plot_temporal_data (gazed_atoms, "All atoms", cfg_iRT_sessionID, cfg_iRT_taskID, axis_sizes=0, ytick_val=0, fig_num=[3,1,3]);
  hold on;
  S  = axes( 'visible', 'off', 'title', ["Atoms gazed, gaze radius:",num2str(cfg_gaze_radius),"px - ",num2str(cfg_iRT_sessionID)," ",cfg_iRT_taskID] );
  hold off;

  gazed_atoms_cumulative_tgt = seriesSum(gazed_atoms_Tgt);
  gazed_atoms_cumulative_int = seriesSum(gazed_atoms_Int);
  gazed_atoms_cumulative = seriesSum(gazed_atoms);
  figure;
  plot_temporal_data (gazed_atoms_cumulative_tgt, "Target atoms", cfg_iRT_sessionID, cfg_iRT_taskID, axis_sizes=0, ytick_val=0, fig_num=[3,1,1]);
  plot_temporal_data (gazed_atoms_cumulative_int, "Interactive atoms", cfg_iRT_sessionID, cfg_iRT_taskID, axis_sizes=0, ytick_val=0, fig_num=[3,1,2]);
  plot_temporal_data (gazed_atoms_cumulative, "All atoms", cfg_iRT_sessionID, cfg_iRT_taskID, axis_sizes=0, ytick_val=0, fig_num=[3,1,3]);
  hold on;
  S  = axes( 'visible', 'off', 'title', ["Cumulative atoms gazed, gaze radius:",num2str(cfg_gaze_radius),"px - ",num2str(cfg_iRT_sessionID)," ",cfg_iRT_taskID] );
  hold off;

% sum of gaze time per atom (gaze/atom)
gaze_per_atom_sum_int = gaze_per_atom_sum_tgt = zeros(atom_count,1);
for a=1:atom_count
  gaze_per_atom_sum_int(a) = sum ( gaze_atomInt_dist(:,a) <=cfg_gaze_radius);
  gaze_per_atom_sum_tgt(a) = sum ( gaze_atomTgt_dist(:,a) <=cfg_gaze_radius);
endfor

  figure();
  plot_data ([1:atom_count],"Atoms", gaze_per_atom_sum_tgt *cfg_t, "Target gaze duration", axis_size=0, ytick_val=0, fig_num=[2,1,1])
  plot_data ([1:atom_count],"Atoms", gaze_per_atom_sum_int *cfg_t, "Interactive gaze duration", axis_size=0, ytick_val=0, fig_num=[2,1,2])
    hold on;
  S  = axes( 'visible', 'off', 'title', ["Sum of gaze time per atom, gaze radius:",num2str(cfg_gaze_radius),"px - ",num2str(cfg_iRT_sessionID)," ",cfg_iRT_taskID] );
  hold off;

  % moving average of gaze region
%figure;plot_angDisp_colored_bg (angDisp, cfg_iRT_sessionID, cfg_iRT_taskID, gaze_region, [2,1]);
%plot_gazeRegion_colored_bg_mw_simple (gaze_region_mw, cfg_iRT_sessionID, cfg_iRT_taskID, [2,2]);

% cumulative of gaze region, 3 plots: gaze region sum, ratio, and ratio discounting the outside region.
%[gaze_region_sum,gaze_region_percent ]= get_gaze_region_sum (gaze_region, cfg_t,cfg_iRT_sessionID, cfg_iRT_taskID);
%plot_data_OLD (gaze_region_sum(:,2),"Gaze Int", gaze_region_sum(:,1),"Gaze Ref", cfg_iRT_sessionID, cfg_iRT_taskID,0,0,[2,3,6,20]);
%plot_temporal_yy(gaze_region_sum,"Regional gaze (s)", seriesSum(angDisp),"AD Sum (deg)", cfg_iRT_sessionID, cfg_iRT_taskID, [2,3,3,21]);
%legend("Tgt gaze", "Int gaze", "Outside gaze", "AD Sum.");

% Compute gaze movement (linear, horizontal and vertical only)
if cfg_gaze_movement == true
  [gaze_px_movement,gaze_px_dx,gaze_px_dy] = compute_gaze_px_movement (gaze_px);
  if cfg_plot_gaze_speed == true
    figure(25);
    plot_temporal_data (gaze_px_dx, "Saccade horizontal speed", cfg_iRT_sessionID, cfg_iRT_taskID, 0,0, [4,1,1]);
    plot_temporal_data (gaze_px_dy, "Saccade vertical speed", cfg_iRT_sessionID, cfg_iRT_taskID,  0,0, [4,1,2]);
    plot_temporal_data (gaze_px_movement, "Saccade speed (absolute)", cfg_iRT_sessionID, cfg_iRT_taskID, 0,0, [4,1,3]);
    plot_temporal_data_colored_bg (angDisp, "Angular disparity", cfg_iRT_sessionID, cfg_iRT_taskID, gaze_region, [4,1,4]);
  endif
  % plot x and y coordinates of gaze in time
  if cfg_plot_gaze_movement == true
    figure;
    plot_temporal_data (gaze_px(:,1), "Gaze x-position", cfg_iRT_sessionID, cfg_iRT_taskID, 0,0, [2,1,1]);
    plot_temporal_data (gaze_px(:,2), "Gaze y-position", cfg_iRT_sessionID, cfg_iRT_taskID,  0,0, [2,1,2]);
  endif
  % plot cumulative gaze movement
  if cfg_plot_gaze_movement_cumulative == true
    cumulative_gaze_px_movement = seriesSum(gaze_px_movement);
    %plot_temporal_data (cumulative_gaze_px_movement, "Saccade", cfg_iRT_sessionID, cfg_iRT_taskID, 0,0, [2,3,1,26]);
  endif
endif

% Find nearest atom from gaze.
% nearest_vertice = [Tgt model atom, Int model atom, dist].
% X model atom = 0 if gaze is outside of window
%LOOK FOR cfg_gaze_dist_matrix AND INSERT CODE THERE (TODO)
frame_count = length(gaze_region);
nearest_vertice = zeros(frame_count,3);
for (i=1:frame_count)
  switch (gaze_region(i))
    case {0}  %when gazing outside
      nearest_vertice(i,:) = [0,0,0];
    case {1}  %when gazing Left (Target) model
      [dist,atom_i] = min ( gaze_atomTgt_dist(i,:) );
      nearest_vertice(i,:) = [atom_i,0,dist];
    case {2}  %when gazing Right (Interactive) model
      [dist,atom_i] = min ( gaze_atomInt_dist(i,:) );
      nearest_vertice(i,:) = [0,atom_i,dist];
  endswitch
endfor
nearest_vertice_B = zeros(frame_count,3);
for (i=1:frame_count)
  switch (gaze_region(i))
    case {0}  %outside
      nearest_vertice_B(i,:) = [0,0,0];
    case {1}  %Left (Target) model
      [dist,atom_i] = min ( gaze_atomTgt_dist(i,:) );
      nearest_vertice_B(i,:) = [dist,atom_i,1];
    case {2}  %Right (Interactive) model
      [dist,atom_i] = min ( gaze_atomInt_dist(i,:) );
      nearest_vertice_B(i,:) = [dist,atom_i,2];
  endswitch
endfor

%{
%run the 2 lines below
 gaze_pupil_diameter_mean = task_data(:,13)+task_data(:,14)/2;
 plot_temporal_data (gaze_region_sum, "Cumulative regional gaze", cfg_iRT_sessionID, cfg_iRT_taskID, 0);
 legend("Target model", "Interactive Model", "Outside");

 swp_frame=330; %frame where exploration ends and adjusts begin

%swap_frame: me-[ 7.8, 13.1, 12.1] she-[ 33.0, 24.7 , 35.6 ]
 scatterplot_2x2 (swp_frame, gazemap_center_dispersion_tgt,gazemap_center_dispersion_int, gazemap_center_speed_tgt,gazemap_center_speed_int, "Gaze dispersion","Center speed", angDisp)
%colormaps: time / angDisp / pupilDiameter
 scatterplot_2x2 (swp_frame, gazemap_center_dispersion_tgt,gazemap_center_dispersion_int, gazemap_centers_delta,gazemap_centers_delta, "Gaze dispersion","Int-Tgt Centers distance")
 scatterplot_2x2 (swp_frame, gazemap_center_dispersion_tgt,gazemap_center_dispersion_int, gazemap_centers_delta,gazemap_centers_delta, "Gaze dispersion","Int-Tgt Centers distance", angDisp)
 scatterplot_2x2 (swp_frame, gazemap_center_dispersion_tgt,gazemap_center_dispersion_int, gazemap_centers_delta,gazemap_centers_delta, "Gaze dispersion","Int-Tgt Centers distance",  (task_data(:,13)+task_data(:,14))/2)
%distancia entre centros vs dilatacao pupilar. Cor por qq
 scatterplot_2x2 (swp_frame, gaze_pupil_diameter_mean,gaze_pupil_diameter_mean, gazemap_centers_delta,gazemap_centers_delta, "Pupil diameter","Int-Tgt Centers distance", angDisp)
 colormap "gray"
 scatterplot_2x2 (swp_frame, gaze_pupil_diameter_mean,gaze_pupil_diameter_mean, gazemap_center_dispersion_tgt,gazemap_center_dispersion_int, "Pupil diameter","Gaze dispersion", angDisp)
 colormap "gray"
% colorbar adds the "color legend"
 scatterplot_2x3 (swp_frame, angDisp, angDisp_discrete, gaze_pupil_diameter_mean,gaze_pupil_diameter_mean, "Angular Disparity", "Pupil diameter")
 scatterplot_2x3 (swp_frame, angDisp, angDisp_discrete, gazemap_centers_delta,gazemap_centers_delta, "Angular Disparity", "Int-Tgt Centers distance")
 scatterplot_2x3 (swp_frame, angDisp, angDisp, gazemap_center_dispersion_tgt,gazemap_center_dispersion_int, "Angular Disparity", "Gaze dispersion")
 scatterplot_2x3 (swp_frame, angDisp_discrete, angDisp_discrete, gazemap_center_dispersion_tgt,gazemap_center_dispersion_int, "Angular Disparity (discrete)", "Gaze dispersion")
% angDisp in log
 scatterplot_1x3 (swp_frame, log(angDisp), gazemap_centers_delta, "log(Angular Disparity)", "Int-Tgt Centers distance")
% pupil diameter x gaze dispersion (existe essa corelacao no objeto interativo?)
 scatterplot_2x3 (swp_frame, gazemap_center_dispersion_tgt, gazemap_center_dispersion_int, gaze_pupil_diameter_mean,gaze_pupil_diameter_mean, "Gaze dispersion", "Pupil Diameter (mm)")
 %}
%DONE
endtext = ["Process complete! Merged data from '", cfg_iRT_input_filename, "' and '", cfg_eyeT_input_filename,"'."];
if (cfg_write_merge_output == true)
  endtext = [endtext, "\nMerged data was stored inside '", writeOutput_merge_fname, "'."];
endif
helpdlg (endtext); %help popup about the process completion
