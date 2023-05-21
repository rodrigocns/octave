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

clear -exclusive *_data;

% SETTINGS
% read eyeTracking .xlsx input file
cfg_eyeT_input = false; %slow, aprox. 1 min per 10k lines and 9 columns
cfg_eyeT_input_filename = "raw_eyeT_r.xlsx"; %name of the input file (.xlsx, numbers only, no commas for decimals)
% (only if necessary) append a column of epoch unix data to the eyeT data array
cfg_add_epoch_from_milliseconds = true;
cfg_epoch_ms_col = 3; %time column in milliseconds used to create epoch column
cfg_epoch_anchor = 1682707674981-5656; %estimated unix epoch at time 0 inside eyetracking data.
% read iRT .xlsx input file
cfg_iRT_input = false; %slow, aprox. 25 seconds for 4.5k lines and 8 columns
cfg_iRT_input_filename = "iRT_data.xlsx"; %name of the iRT sheets file unpackaged by unpacking_sheets.m
% merge eyeTracker data to the chosen iRT task (be sure that the files are from the same session)
cfg_data_merge = true;
cfg_iRT_sessionID = 1682707472090; %session ID of the desired task
cfg_iRT_taskID = "bolaBastao_c"; %task ID of the desired task
cfg_iRT_cols = [3:8]; %range of desired data columns from raw_iRT_data. 5:8 is quaternion data, 3 is unix epoch
cfg_eyeT_cols = [1,2,4,7:10]; %range of desired data columns from eyeT_data. epoch data was appended as the first column
% write output file from task_data
cfg_write_output = true;
cfg_output_filename = strcat("mergeOutput_", cfg_iRT_taskID, num2str(cfg_iRT_sessionID), ".xlsx"); % name of the output file (.xlsx) from the merging of eyeT_data and iRT_data

%FUNCTIONS
% recognize input file format and reads eyeTracker data
function [raw_eyeT_data, eyeT_header_data] = eyeT2oct (filename)
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
%add epoch data to a time column counting milliseconds
function epoch_column = add_epoch (raw_eyeT_data, epoch_column, epoch_anchor)
  printf("Adding epoch column..");
  epoch_column = raw_eyeT_data(:,epoch_column);
  epoch_column = epoch_column + epoch_anchor;
  printf("Done!\n");
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
%==========================

%SCRIPTS
%data checks!
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
  [raw_eyeT_data, eyeT_header_data] = eyeT2oct (cfg_eyeT_input_filename);
else
  disp("Skipping eyeT file read.");
endif
%eyeTracking data pre-process: add epoch column to eyeT_data and eyeT_header_data
if cfg_add_epoch_from_milliseconds == true
  eyeT_data = horzcat( add_epoch (raw_eyeT_data, cfg_epoch_ms_col, cfg_epoch_anchor) , raw_eyeT_data );
  eyeT_header_data = horzcat ('EyeT Epoch', raw_eyeT_header_data);
endif
%iRT data input
if cfg_iRT_input == true
  [session_data,raw_iRT_data] = iRT2oct (cfg_iRT_input_filename);
else
  disp("Skipping iRT file read.");
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
% write output file
if cfg_write_output == true
  writeOutput (cfg_output_filename, merged_header_data, task_data);
endif
%clear cfg variables for cleaner debugging
clear cfg*
