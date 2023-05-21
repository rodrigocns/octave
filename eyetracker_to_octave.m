pkg load io
%{
This script atempts to read the eyeTracker and iRT files to pre-process them
and then merge both data in time. We used the UNIX EPOCH time to do this, as it
is the innate time used in most programming languages, especially javascript and
HTML. If the data from your eyetracker device/software does not contain this
kind of information, then it must be added using the epoch data produced from
iRT, or changing from another epoch system into the unix epoch.
%}

clear -x raw* session_data iRT_data

% SETTINGS
% read eyeTracking .xlsx input file
cfg_eyeT_input = false;
cfg_eyeT_input_filename = "raw_eyeT_r.xlsx"; %name of the input file (.xlsx, numbers only, no commas for decimals)
% append a column of epoch unix data to the eyeT data array (only if necessary)
cfg_add_epoch_from_milliseconds = true;
cfg_epoch_ms_col = 3;%time column in milliseconds used to create epoch column
cfg_epoch_anchor = 1682707674981-5656; %estimated unix epoch at time 0.
% read iRT .xlsx input file
cfg_iRT_input = false;
cfg_iRT_input_filename = "iRT_data.xlsx"; %name of the iRT sheets file unpackaged by unpacking_sheets.m
% merge eyeTracker data to the chosen iRT task (be sure that the files are from the same session)
cfg_iRT_merge = true;
cfg_iRT_sessionID = 1682707472090; %session ID of the desired task
cfg_iRT_taskID = "bolaBastao_c"; %task ID of the desired task
cfg_iRT_cols = 3:8; %range of desired data columns from task data. 5:8 is quaternion data, 3 is unix epoch

%FUNCTIONS
% recognize input file format and reads eyeTracker data
function raw_eyeT_data = eyeT2oct (filename)
  tic();
  printf("Reading eyetracker data (this can take a while)..");
  raw_eyeT_data = xlsread (filename);
  printf(".done!");
  toc();
endfunction

%add epoch data to a time column counting milliseconds
function epoch_column = add_epoch (raw_eyeT_data, epoch_column, epoch_anchor)
  printf("Adding epoch column..");
  epoch_column = raw_eyeT_data(:,epoch_column);
  epoch_column = epoch_column + epoch_anchor;
  printf("Done!\n");
endfunction
% read data from iRT .xlsx. session_arr has session details, raw_arr has the data bulk in a cell matrix
function [session_arr,raw_arr] = iRT2oct (filename)
  tic();
  printf("Reading iRT data (this can take a while)..");
  xls = xlsopen(filename);
  session_arr = xls2oct(xls,"sessions");
  raw_arr = xls2oct(xls,"data");
  [xls] = xlsclose(xls);
  printf(".done!");
  toc();
endfunction
% create data array and its header from slice of iRT data (task data) corresponding to session_ID and task_ID.
% 'range_of_columns' is the range of the desired columns from raw_iRT_data (5:8 or 3,4,5:8)
function [matrix_quat_data,data_header] = slice_task_data (raw_arr, session_ID, task_ID, range_of_columns)
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
  matrix_quat_data = cell2mat( raw_arr( first_line:last_line, range_of_columns) );
  data_header = cell(raw_arr(1,range_of_columns));
  printf(".data sliced!\n");
endfunction
%==========================
%SCRIPTS
%eyetracking data input
if cfg_eyeT_input == true
  raw_eyeT_data = eyeT2oct (cfg_eyeT_input_filename);
else
  disp("Skipping eyeT file read. \n");
endif

%eyeTracking data pre-process
if cfg_add_epoch_from_milliseconds == true
  eyeT_data = horzcat( raw_eyeT_data, add_epoch (raw_eyeT_data, cfg_epoch_ms_col, cfg_epoch_anchor) );
endif

%iRT data input
if cfg_iRT_input == true
  [session_data,raw_iRT_data] = iRT2oct (cfg_iRT_input_filename);
else
  disp("Skipping iRT file read. \n");
endif

%iRT - eyeTracking data merge
if cfg_iRT_merge == true
  task_data = slice_task_data (raw_iRT_data, cfg_iRT_sessionID, cfg_iRT_taskID, cfg_iRT_cols);

endif

%clear cfg variables
clear cfg*
