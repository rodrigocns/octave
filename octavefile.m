pkg load io
#{ 
#json method
json_string = fileread('OctaveAccessTest.json');
data_obj = fromJSON(json_string);
#}
tic;
sheets_filename = 'Octave32767str.xlsx';
toc;
[sheets_num,sheets_txt,sheets_raw] = xlsread(sheets_filename,1);
split(sheets_raw[],",");