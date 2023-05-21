pkg load io

clear -x raw_eyeT_data
% settings

config_tipo = 1; %1=cor || 2=cinza || qual o tipo do teste.
config_dados_input = 1; %0=nao || 1=sim || importar dados do xlsx? so precisa 1 vez.
config_table_output = 1; %0=nao || 1=sim || imprimir arquivo xlsx dos dados?
config_plot = 1; %0=nao || 1=sim || plotar os dados da tabela


function time_calculated = calc_t(time_iRT,trio_valores)
  %%converte tempo distorcido do iRoT para tempo correto do opengaze com base no
  %% tempo decorrido/atual do iRot e do tempo inicial e duracao do opengaze
  dur_iRT      = trio_valores(1);
  t0_eyeT       = trio_valores(2);
  dur_eyeT  = trio_valores(3);
  time_calculated = (time_iRT*dur_eyeT/dur_iRT)+t0_eyeT;
  %formulae from excel (eyeT=eyeTracker, dur=duration, iRT=interactiveRotationTest):
  % time_calculated = (t - t0_iRoT) * (dur_eyeT/duration_iRT) + t0_et
  % time_calculated = (t - 71) * (86.8065/79.9) + 40.8895
  %as t is already subtracted from t0_iRot,
  % time_calculated = (0) * (86.8065/79.9) + 40.8895
endfunction
function linha_de_dados = calcular_linha(id_linha,moving_window)
  i=5; r=2;
  %{
   'moving_window' is a sample of data from used devices like eye-tracker,
   expected to have higher data frequencies (such as 60Hz) than iRT (10Hz).
   'i' is the index of the data sample. It is a sample of 9 points going from
   i-4 to i+4. Internally, octave will see 'i' going from 1 to 9, as
   i=5 is the middle of the sample array, expected to be the closest in time to
   the iRT data. 'r' is the range used to calculate means or moving averages.

   Each column of eye-tracker data could hold distinct types of data, and each
   one of these types should be interpolated acordingly.
   Here are some suggestions (evaluate each one of your data columns):

   Continuous real data: mean of i (from -r to +r)
   Ex.: X coordinates of the eye position as a % of the screen
   Continuous integer data: round the mean of i (from -r to +r)
   Ex: epoch/timetick/system time count
   id or itemized: take i itself
   Ex: eye fixation count
   reference : take i itself
   Ex: initial time of fixation POG
   event : take the event happening (max value from entire window)
   Ex: mouse state (when it was pressed or released)
  %}
  linha_de_dados(1,1) = moving_window(i,1); %CNT
  linha_de_dados(1,2) = mean (moving_window(i-r:i+r,2) ); %time
  linha_de_dados(1,3) = round (mean (moving_window(i-r:i+r,3) ) ); %timetick
  linha_de_dados(1,4) = mean (moving_window(i-r:i+r,4) ); %FPOGX
  linha_de_dados(1,5) = mean (moving_window(i-r:i+r,5) ); %FPOGY
  linha_de_dados(1,6) = moving_window(i,6); %FPOGS
  linha_de_dados(1,7) = moving_window(i,7); %FPOGD
  linha_de_dados(1,8) = moving_window(i,8); %FPOGID
  linha_de_dados(1,9) = moving_window(i,9); %FPOGV
  linha_de_dados(1,10) = mean (moving_window(i-r:i+r,10) ); %BPOGX
  linha_de_dados(1,11) = mean (moving_window(i-r:i+r,11) ); %BPOGY
  linha_de_dados(1,12) = moving_window(i,12); %BPOGV
  linha_de_dados(1,13) = mean (moving_window(i-r:i+r,13) ); %CX
  linha_de_dados(1,14) = mean (moving_window(i-r:i+r,14) ); %CY
  linha_de_dados(1,15) = max (moving_window(:,15) ); %CS
  linha_de_dados(1,16) = mean (moving_window(i-r:i+r,16) ); %LPCX
  linha_de_dados(1,17) = mean (moving_window(i-r:i+r,17) ); %LPCY
  linha_de_dados(1,18) = mean (moving_window(i-r:i+r,18) ); %LPD
  linha_de_dados(1,19) = mean (moving_window(i-r:i+r,19) ); %LPS
  linha_de_dados(1,20) = moving_window(i,20); %LPV
  linha_de_dados(1,21) = mean (moving_window(i-r:i+r,21) ); %RPCX
  linha_de_dados(1,22) = mean (moving_window(i-r:i+r,22) ); %RPCY
  linha_de_dados(1,23) = mean (moving_window(i-r:i+r,23) ); %RPD
  linha_de_dados(1,24) = mean (moving_window(i-r:i+r,24) ); %RPS
  linha_de_dados(1,25) = moving_window(i,25); %RPV
  linha_de_dados(1,26) = moving_window(i,26); %BKID
  linha_de_dados(1,27) = max (moving_window(:,27) ); %BKDUR
  linha_de_dados(1,28) = moving_window(i,28); %BKPMIN
  linha_de_dados(1,29) = mean (moving_window(i-r:i+r,29) ); %LPMM
  linha_de_dados(1,30) = moving_window(i,30); %LPMMV
  linha_de_dados(1,31) = mean (moving_window(i-r:i+r,31) ); %RPMM
  linha_de_dados(1,32) = moving_window(i,32); %RPMMV
  linha_de_dados(1,33) = max (moving_window(:,33) ); %SACCADE_MAG
  linha_de_dados(1,34) = max (moving_window(:,34) ); %SACCADE_DIR
endfunction



%%inicializar valores
%valores_ref = [duracao iRoT, t0 opengaze, duracao opengaze]
val_ref(1,1:3) = [79.8, 40.8895, 86.8065]; %cor
val_ref(2,1:3) = [81.2, 131.384, 88.506]; %cinza
if (config_tipo == 1)
  nome_aba_output = "cor";
elseif (config_tipo == 2)
  nome_aba_output = "cinza";
else
  error("ERRO! config_tipo nao reconhecido!");
endif
%%inicializa tabela(n x m) de dados com n linhas de zeros
% n=linhas de dados do iRT (duracao*10Hz+1)
% m = categorias de dados do eyetracker
tabela(1:( val_ref(config_tipo,1) * 10 + 1 ) , 1:35) = 0;
% preenche ultima coluna com dados de tempo do iRT (tempo avancando a cada 0.1s)
tabela(1:end,35) = [0:0.1:val_ref(config_tipo,1)];

%%abre dados da tabela openGaze
if (config_dados_input == 1)
  tic();
  printf ("Carregando dados (pode demorar)...");
  raw_eyeT_data = xlsread ("User 1_all_gaze.xlsx", 1, "C3:AJ13843");
  %PODE SER TUDO SE FOR CELL e depois a pessoa escolhe quais colunas interessam
  printf ("pronto!");
  toc();
endif

%merging loop script
tic();
for i = 1 : size(tabela,1)
  %achar posicao do melhor valor de tempo
  [x,dados_count] = min( abs( raw_eyeT_data(:,2)-calc_t (tabela(i,35),val_ref(config_tipo,1:3)) ) );
  %%preencher tabela com valores da linha encontrada
  tabela(i,1:34) = calcular_linha (dados_count, raw_eyeT_data(dados_count-4:dados_count+4,:));
  if (rem(i,100) == 0) %remainder of division
    printf("[%i]",i);
  endif
endfor
toc();

% Write .xlsx output file
tic();
if (config_table_output == 1)
  printf(cstrcat("Gravando dados da tabela ",nome_aba_output," ..."));
  xlswrite ("tabela_de_dados_opengaze.xlsx", tabela, nome_aba_output, "A2");
  printf("feito!");
endif
toc();

% plot (?)
if (config_plot==1)
  i=2;
  while (i <= size(tabela,1))
    plot3 (tabela(i-1:i,4), tabela(i-1:i,5), tabela(i-1:i,2));
  %  plot3 (tabela(1:100,4), tabela(1:100,5), tabela(1:100,2));
    hold on;
    xlabel("x"); ylabel("y"); zlabel("time");
    axis ([-1 1 -1 1 0 tabela(end,2)])
    pause(0.1) %tempo para cada passo
    i++;
  endwhile
endif

%%EOF
