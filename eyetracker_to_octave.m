pkg load io
%{
This script atempts to read the eyeTracker and iRT files to pre-process them
and then merge both data in time. We used the UNIX EPOCH time to do this, as it
is the innate time used in most programming languages, especially javascript and
HTML. If the data from your eyetracker device/software does not contain this
kind of information, then it must be added using the epoch data produced from
iRT, or changing from another epoch system into the unix epoch.
%}


clear -x raw_eyeT_data
% SETTINGS
% read .xlsx input file
cfg_dados_input = false;
filename_raw_eyeT = "raw_eyeT_r.xlsx"; %name of the input file (.xlsx, numbers only, no commas for decimals)
% add epoch data to a time column counting milliseconds (only if necessary)
cfg_add_epoch_from_milliseconds = true;
cfg_epoch_ms_col = 3;%time column in milliseconds used to create epoch column
cfg_epoch_anchor = 1682707674981-5656; %estimated unix epoch at time 0.
%

%FUNCTIONS
% recognize input file format and reads eyeTracker data
function raw_eyeT_data = eyeT2oct (filename)
  tic();
  printf("Reading .xlsx data..");
  disp("(this can take a while)");
  raw_eyeT_data = xlsread (filename);
  disp(".done!");
  toc();
endfunction

%add epoch data to a time column counting milliseconds
function epoch_column = add_epoch (raw_eyeT_data, epoch_column, epoch_anchor)
  disp("Adding epoch column..");
  epoch_column = raw_eyeT_data(:,epoch_column);
  epoch_column = epoch_column + epoch_anchor;
  disp("Done!\n");
endfunction
%SCRIPTS
%eyetracking data input
if cfg_dados_input == true
  raw_eyeT_data = eyeT2oct (filename_raw_eyeT);
else
  disp("Skipping file read. \n");
endif

%eyeTracking data pre-process
if cfg_add_epoch_from_milliseconds == true
  eyeT_data = horzcat( raw_eyeT_data, add_epoch (raw_eyeT_data, cfg_epoch_ms_col, cfg_epoch_anchor) );
endif

%iRT data input

