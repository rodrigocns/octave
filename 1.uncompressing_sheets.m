pkg load io
%Script lê arquivo 'input_filename' e cria um arquivo csv 'output_filename'
% expandindo de uma resolução por i_line para um ponto de dado por i_line.

clear
% Settings (altere conforme necessário) ====================
input_filename = 'iRT gsheets.xlsx'; %downloaded sheets filename (format incl.).
output_filename = 'iRT_data'; %sem o formato de arquivo
zipped_rows = [8,9,10,11,12,13]; %colunas do input para descomprimir
header={"sessionID", "task_id", "epoch", "duration", "Qi", "Qj", "Qk", "Qr"};
header_atemporal={"sessionID", "task_id", "startEpoch", "email", "subject", "pxAngsRatio"};
% Funções =================================================
function uncompressed_data = extract(compressed_data)
  splitted_data = strsplit(compressed_data, ',');
  uncompressed_data = str2double(splitted_data)';
endfunction

% Copiar dados da planilha para o octave =================
tic;printf("Lendo planilha...");
[sheets_num,sheets_txt,sheets_raw] = xlsread(input_filename,1);
printf("DONE! ");
toc;

% Descompactando dados ===================================
tic;printf("Descomprimindo dados...");

full_matrix = header;
atemporal_cmatrix = header_atemporal;
%para cada linha de dados compactos (2a até a ultima)
contador = 0;
for i_line = 2:size(sheets_raw,1)
  % pega dados únicos da sessão/task (session_id, task_id, etc)
  session_id = sheets_raw{i_line,4};
  task_id = sheets_raw{i_line,5};
  pxAngsRatio = mean(extract(sheets_raw{i_line,6}));

  % descompactar cada array de dados na matriz temporal_array (de uma para várias linhas)
  % uma coluna por vez
  temporal_array = [];
  for i_col = 1:size(zipped_rows,2)
    temporal_array = [temporal_array, extract( sheets_raw{ i_line, zipped_rows(i_col) }) ];
  endfor

  % epoch_n recebe epoch_0 + delta (linha i_line, coluna 7)
  temporal_array(:,1) = temporal_array(:,1) + str2num(sheets_raw{i_line,7});

  % cria cell matrix n x 2, n -> comprimento de temporal_array
  task_cmatrix = cell(size(temporal_array,1),2);
  task_cmatrix(:,1) = {session_id};  % preenche IDs
  task_cmatrix(:,2) = {task_id};
  task_cmatrix = horzcat(task_cmatrix,num2cell(temporal_array));

  % append na matriz completa
  full_matrix = vertcat(full_matrix,task_cmatrix);

  % append na matriz de dados constantes da sessao
  atemporal_cmatrix = vertcat(atemporal_cmatrix, {session_id, task_id, sheets_raw{i_line,7}, sheets_raw{i_line,2}, sheets_raw{i_line,3}, pxAngsRatio });

  contador++;
endfor
printf("%i linhas processadas. DONE! ",contador);toc;
clear contador;

% Salvando os dados ====================================
tic;printf("Salvando dados...");

% write to xlsx file
xls = xlsopen(strcat(output_filename,".xlsx"), 1);
xls = oct2xls(full_matrix,xls,"data");
xls = oct2xls(atemporal_cmatrix,xls,"sessions");
xls = xlsclose(xls);

printf(".xlsx salvo. ");

printf("DONE! ");toc;

% CSV============================

% Write to CSV file
%cell2csv(strcat(output_filename,".csv"), full_matrix);
%printf(".csv salvo ");

#{
tic;printf("Carregando dados...");
cell_mat_new = csv2cell(strcat(output_filename,".csv"));

printf("DONE! ");toc;
#}

%{
fid = fopen('iRT.bin', 'r');
str_data = fread(fid, [n, 2], 'string');
num_data = fread(fid, [n, 6], 'double');
fclose(fid);
%}
