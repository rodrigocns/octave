pkg load io
%Script lê arquivo 'input_filename' e cria um arquivo csv 'output_filename'
% expandindo de uma resolução por linha para um ponto de dado por linha.

clear
input_filename = 'iRT gsheets.xlsx'; %downloaded sheets filename (format incl.).
output_filename = 'iRT_data'; %sem o formato de arquivo
colunas = [8,9,10,11,12,13]; %colunas do input para descomprimir

% Funções =================
function data_ready = extract(compressed_data)
  splitted_data = strsplit(compressed_data, ',');
  data_ready = str2double(splitted_data)';
endfunction

% Copiar dados da planilha ============================
tic;printf("Lendo planilha...");
[sheets_num,sheets_txt,sheets_raw] = xlsread(input_filename,1);
printf("DONE! ");
toc;

% Descompactando dados ===================================
tic;printf("Descomprimindo dados...");
header={"sessionID", "task_id", "epoch", "duration", "Qi", "Qj", "Qk", "Qr"};
full_matrix = header;

%para cada linha de dados (2a até a ultima)
contador = 0;
for linha = 2:size(sheets_raw,1)
  numeric_matrix = [];
  %para cada coluna de dados compact.
  for col_id = 1:size(colunas,2);
    %descompactar dados das celulas em uma matriz inteira
    numeric_matrix = [numeric_matrix,extract(sheets_raw{linha,colunas(col_id)})];
  endfor
  %soma epoch inicial na coluna de epoch
  numeric_matrix(:,1) = numeric_matrix(:,1) + str2num(sheets_raw{linha,7});
  %criando uma cell matrix para os dados da task
  session_id = sheets_raw{linha,4};
  task_id = sheets_raw{linha,5};
  cell_matrix = cell(size(numeric_matrix,1),2);  %cria cell matrix n x 2
  cell_matrix(:,1) = {session_id};  %preenche IDs
  cell_matrix(:,2) = {task_id};
  task_matrix = horzcat(cell_matrix,num2cell(numeric_matrix));
  full_matrix = vertcat(full_matrix,task_matrix);
  contador++;
endfor
printf("%i linhas processadas. DONE! ",contador);toc;
clear contador;

% Salvando os dados ====================================
tic;printf("Salvando dados...");



% Write to XLSX file (single line)
%oct2xls(strcat(output_filename,".xlsx"), full_matrix);
%xlswrite (strcat(output_filename,".xlsx"), full_matrix);

% write to xlsx file
xls = xlsopen(strcat(output_filename,".xlsx"), 1);
xls = oct2xls(full_matrix,xls,1);
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
