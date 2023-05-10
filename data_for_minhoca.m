pkg load io
clear;
% Find data for computing the "quaternionic worm-path"
%from the "filename" excel spreadsheet file of unpackaged iRT data from gsheets


session_ID = 1682699553789;
task_ID = "bolaBastao_c";
filename = "iRT_data.xlsx";

% read data from .xlsx
xls = xlsopen(filename);
session_arr = xls2oct(xls,"sessions");
raw_arr = xls2oct(xls,"data");
[xls] = xlsclose(xls);

printf("Reading done\n");

% get reference data

for n = 1 :size(session_arr,1)
  if and( strcmp(session_arr{n,1},num2str(session_ID)), strcmp(session_arr{n,2}, task_ID) )
    ref_quat_data = cell2mat(session_arr(n,7:10));
    printf("Reference orientation found: %d %d %d %d\n",ref_quat_data(1),ref_quat_data(2),ref_quat_data(3),ref_quat_data(4));
    break;
  endif
endfor

% find slice of data corresponding to session_ID and task_ID
% find first line of the slice
first_line= -1;
for n = 1 :size(raw_arr,1)
  if and( strcmp(raw_arr{n,1},num2str(session_ID)), strcmp(raw_arr{n,2}, task_ID) )
    first_line = n;
    printf("First line of slice is %i\n", first_line);
    break;
  endif
endfor
% find last line of the slice
last_line= -1;
for n = first_line : size(raw_arr,1)
  if not(and( strcmp(raw_arr{n,1},num2str(session_ID)), strcmp(raw_arr{n,2}, task_ID) ))
    last_line = n-1;
    printf("Last line of slice is %i\n", last_line);
    break;
  endif
endfor

% get slice of data
matrix_quat_data=cell2mat(raw_arr(first_line:last_line,5:8));
clear -x matrix_quat_data ref_quat_data;
%quaternion : [i j k r]
printf("DONE!\n");

