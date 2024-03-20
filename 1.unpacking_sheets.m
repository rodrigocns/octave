pkg load io
%{
 Script reads file 'input_filename' containing the packaged data downloaded from
 gsheets, unpackage the session data contained in each line to the number of
 lines it should have (1 line per 0.1 second), and writes the unpackaged data in
 a .xlsx file named 'output_filename'.
%}
clear

% Settings (modify as needed) ====================

% boolean about the file origin: is the input file a local backup data file or a packaged data file?
is_local_backup = false;

% name of downloaded sheets filename (default)
input_filename = 'iRT gsheets.xlsx';

% output file name (default)
output_filename = 'iRT data.xlsx';

% Header names of the data tab (temporal data)
header={"sessionID", "task_id", "epoch", "duration", "Qi", "Qj", "Qk", "Qr"};
% Header names of the session identifiers tab (ex.: name, id, email, etc.)
header_atemporal={"sessionID", "task_id", "startEpoch", "email", "subject", "pxAngsRatio"};
% The columns from input file that should be expanded.
zipped_rows = [8,9,10,11,12,13];

%{
   #=========================================#
   # DON'T MODIFY ANYTHING BELLOW THIS LINE! #
   #=========================================#
%}
% functions =================================================
function unpackaged_data = extract (packaged_data)
  splitted_data = strsplit (packaged_data, ',');
  unpackaged_data = str2double (splitted_data)';
endfunction

% Prompt user for filenames ==================
input_prompt = inputdlg ("Insert the name of the packaged data file previously downloaded. Ex: iRT gsheets.xlsx", "Input data file");
if (isempty ( char (input_prompt) ) == 1 )
  printf ( cstrcat ("Blank input, using default file name: '",input_filename,"'.\n"));
else
  if ( endsWith (char (input_prompt), ".xlsx" ) == 0)
    warndlg ( cstrcat ("File name '",char (input_prompt), "' does not end with .xlsx\nAre you sure it is in the right file format? "), "Input data file");
    endif
  input_filename = char(input_prompt);
endif

output_prompt = inputdlg ("Insert the name of the unpackaged data file to be written. Ex: iRT data.xlsx", "Output data file");
if (isempty ( char (output_prompt) ) == 1 )
  printf ( cstrcat ("Blank input, using default file name: '",output_filename,"'.\n"));
else
  if ( endsWith (char (output_prompt), ".xlsx" ) == 0)
    warndlg ( cstrcat ("File name '",char (output_prompt), "' does not end with .xlsx\nAre you sure it is in the right file format? "), "Input data file");
    endif
  output_filename = char(output_prompt);
endif

if (strcmp (input_filename, output_filename) == 1 )
  output_filename = strcat( "unpacked ",output_filename);
  warndlg ( cstrcat ("Input name is the same as output name. Changed output name to '",output_filename, "' to avoid overwrites.\n"));
endif

is_back_prompt = inputdlg ("Is the input file a local backup file?\nType Y for yes, N for no.", "Local backup?");
if (is_back_prompt{1} == "Y" )
  printf(cstrcat ("The file '", input_filename, "' has been declared as a local backup and will be split instead of unpackaged.\n"));
  is_local_backup = true;
elseif (is_back_prompt{1} == "N" )
  printf(cstrcat ("The file '", input_filename, "' has been declared as a packaged file and will be unpackaged.\n"));
  is_local_backup = false;
else
  printf ( cstrcat ("Blank or unrecognized input. The file '", input_filename, "' will be considered as a packaged data file.\n"));
endif

% Read sheets data to octave =====================
tic;printf ("Reading Sheets...");
[sheets_num,sheets_txt,sheets_raw] = xlsread (input_filename,1);
printf("DONE! ");
toc;

% Unpacking data (if input file is not a local backup) ==========
if (is_local_backup == false)
  tic;printf ("Unpacking data...");

  temporal_matrix = header;
  atemporal_cmatrix = header_atemporal;
  %for each line of packaged data (2nd till last)
  counter = 0;
  for i_line = 2:size (sheets_raw,1)
    % take unique data from session/task (session_id, task_id, etc)
    session_id = sheets_raw{i_line,4};
    task_id = sheets_raw{i_line,5};
    pxAngsRatio = mean (extract (sheets_raw {i_line,6} ) );

    % unpackage each data array in temporal_array matrix (from a single line to many lines)
    % one column at a time
    temporal_array = [];
    for i_col = 1:size (zipped_rows,2)
      temporal_array = [temporal_array, extract( sheets_raw{ i_line, zipped_rows(i_col) }) ];
    endfor

    % epoch_n receives epoch_0 + delta (line 'i_line', column 7)
    temporal_array(:,1) = temporal_array(:,1) + str2num(sheets_raw{i_line,7});

    % creates cell matrix n x 2 for sessionID and taskID, n -> length of temporal_array
    task_cmatrix = cell(size(temporal_array,1),2);
    task_cmatrix(:,1) = {session_id};  % fill in IDs
    task_cmatrix(:,2) = {task_id};
    task_cmatrix = horzcat(task_cmatrix,num2cell(temporal_array));

    % append in complete matrix
    temporal_matrix  = vertcat(temporal_matrix ,task_cmatrix);

    % append in matrix of constant session data
    atemporal_cmatrix = vertcat(atemporal_cmatrix, {session_id, task_id, sheets_raw{i_line,7}, sheets_raw{i_line,2}, sheets_raw{i_line,3}, pxAngsRatio });

    counter++;
  endfor
  printf("%i lines processed. DONE! \n",counter);toc;
  clear counter;
endif

% Splitting data (if file is local backup) ==================
if (is_local_backup == true)
  tic;printf ("Splitting data in two sheets...");

  % Extract columns A to H for the "data" sheet
  temporal_matrix  = sheets_raw(:, 1:8);

  % Extract columns I to AD for the "sessions" sheet
  atemporal_cmatrix = sheets_raw(:, 9:30);
  printf(" DONE! \n");toc;
endif

% Saving data ====================================
tic;
printf ( cstrcat ( "Saving ", num2str(size(temporal_matrix ,1)), " rows of data (estimated time: ", num2str(size(temporal_matrix,1)/40), " s)..."));
% write to xlsx file
xls = xlsopen(strcat(output_filename), 1);
xls = oct2xls(temporal_matrix ,xls,"data");
xls = oct2xls(atemporal_cmatrix,xls,"sessions");
xls = xlsclose(xls);

printf(strcat(output_filename," saved. "));

printf("DONE! \n");toc;
helpdlg ( cstrcat ("Unpackaging of file '",input_filename, "' is complete.\nData was stored inside '", output_filename, "'.") );
