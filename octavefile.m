pkg load io
#sheets_filename = 'Octave32767str.xlsx';
sheets_filename = 'iRT gsheets.xlsx';
col_to_decompress = [6,7,8,9,10,11]; 

tic;
[sheets_num,sheets_txt,sheets_raw] = xlsread(sheets_filename,1);
toc;

function data_ready = decompress_data(compressed_data)
  splitted_data = strsplit(compressed_data, ',');
  data_ready = str2double(splitted_data);
endfunction



teste = decompress_data(sheets_raw{2,3});
teste
#{
sheets_raw

sheets_raw{:,1}



#}