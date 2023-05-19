pkg load io
pkg load quaternion
%http://jmol.sourceforge.net/demo/jssample0.html

%{
 v0.9 - get ready for the manuscript
 Dados da tarefa (task) agora estao dentro do arquivo de funcao get_task_data.m
%}
% ===CONFIGS===
task = ["C","G","I"];      %identifiers of each task. Rename accordingly.
config_task =           2; %Refer to line above. Usado para obter dados referentes a task
config_input_data =     1; %1 para ler os .xlsx de input
config_pre_calc =       1; %1 se for realizar outros calculos (importantes para o resto)
config_convhull =       0; %1 para calcular a convexhull da projecao xy dos atomos em cada tempo
config_gaze_atom_dist = 1; %1 para calcular matriz de distancia entre gazepoint e atomos interativos.
config_dist_integral =  1; %1 para calcular coluna de alfa/transparencia de cada atomo menos visto
config_dist_integral_temporal =  1; %1 para calcular coluna de alfa/transparencia de cada atomo menos visto
config_xls_write =  "0"; %1: convexhull, 2:gazepoint. Escrever arquivos.
config_graficos =   "0"; %1:projecao, 2:projecao em t, 3:scatter3,  4:scatter3 em frame, 5:scatter3 em frame com gaze mais proximo
subject =  63; %1~62
config_screen_size = [1360,720];
config_ref_center_px = [211,376];
config_cvs_center_px = [615,379];
config_fator_px = 23.699; %razao entre fracao da tela e posicao em px (calculo em User 1_all_gaze.xlsx > "px to x coord")
config_gauss_wdt = 0.922; %desvio padrao do config_fator_px
%xlsx_file_name = strcat ("compiladoPorSubjects_", task(config_task), ".xlsx");
xlsx_file_name = strcat ("tabela_de_dados_opengaze",".xlsx");
[config_model_name,task_name,config_ref_quat] = get_task_data (config_task);
xyz_file_name = strcat ("modelos/", config_model_name, ".xyz"); %name of .xyz file used in task


% ===FUNCTIONS===

% rot_matrix retorna matriz de rotacao {R} a partir de um quaternio {qr,qi,qj,qk}
function R = rot_matrix (qr,qi,qj,qk, s = 1) %no .xlsx esta como qi,qj,qk,qr.
  R= [1-2*s*(qj^2 + qk^2), 2*s*(qi*qj - qk*qr), 2*s*(qi*qk + qj*qr);
      2*s*(qi*qj + qk*qr), 1-2*s*(qi^2 + qk^2), 2*s*(qj*qk - qi*qr);
      2*s*(qi*qk - qj*qr), 2*s*(qj*qk + qi*qr), 1-2*s*(qi^2 + qj^2)];
endfunction

% axis-angle to quaternions || in:(x,y,z,angle in degrees) || out:[qw,qx,qy,qz] (qw is real part)
function Q = axangle2quat (x,y,z,angle) %Q = [qw,qx,qy,qz].
  Q = [ cos(deg2rad(angle/2)), x*sin(deg2rad(angle/2)), y*sin(deg2rad(angle/2)), z*sin(deg2rad(angle/2))];
endfunction

% codify atoms by element|| in:(atom_count,atom_xyz,atom_elem) || out:[size,R,G,B]
%used in some graph renderings at the end of this script
function atom_cor = generate_color_vector (atom_count, atom_xyz, atom_elem)
  for i = 1:atom_count
    switch ( strvcat(atom_elem(i)) )  %strvcat extracts a string from a cell array
      case "H"
        atom_xyz(i,4) = 1; % atomic number
        atom_cor(i,1:4) = [74, 0.75,0.75,0.75]; %[size, R,G,B]
      case "C"
        atom_xyz(i,4) = 12;
        atom_cor(i,1:4) = [154, 0,0,0];
      case "N"
        atom_xyz(i,4) = 14;
        atom_cor(i,1:4) = [140, 0,0,1];
      case "O"
        atom_xyz(i,4) = 16;
        atom_cor(i,1:4) = [146, 1,0,0];
      case "Co"
        atom_xyz(i,4) = 27;
        atom_cor(i,1:4) = [100, 1,0,0];
      otherwise
        error ("ERRO! Novo elemento detectado. Atualizar codigo");
        atom_xyz(i,4) = -1;
        atom_cor(i,1:4) = [20, 0,1,0];
    endswitch
  endfor
endfunction

% normalize xyz atom coords 0,0,0 (at the center of jmol boundingbox )
function norm_atom_xyz = normalize_jmol_rot_center (atom_xyz)
  % get max and min x,y,z coords from the atom_xyz array.
  % These are the boundingbox extremities in jmol.
  max_xyz = max(atom_xyz(:,1:3));
  min_xyz = min(atom_xyz(:,1:3));
  % normalize position of entire array
  normalization_center = (max_xyz+min_xyz)/2;
  norm_atom_xyz = atom_xyz(:,1:3) - normalization_center;
endfunction

% returns cell array with element symbols and xyz coordinates matrix(nx3) of atoms in .xyz file
function [atom_count,elem,atom_coords] = get_xyz_data (filename)
  printf(strcat("Opening .xyz file: ",filename,"... "));
  fid = fopen(filename, 'r');
  % Read the number of atoms from the first line of the file
  printf(" reading...");
  line = fgetl(fid);
  atom_count = str2num(line);

  % Pre-allocate arrays for the element symbols and coordinates
  elem = cell(atom_count, 1);
  atom_coords = zeros(atom_count, 3);

  % Loop over each line of the file after the first line
  line = fgetl(fid);
  for i = 1:atom_count
    line = fgetl(fid);
    line_data = strsplit(line);
    elem{i} = line_data{1};
    atom_coords(i, :) = str2double(line_data(2:4));
  end
  fclose(fid);
  printf(" sucess! \n");
endfunction



% ===== Data Input =====
% Obtaining data

if (config_input_data == 1)
  %clean all variables except those below
  clear -x config* subject *name
  tic();
  %take count, element and xyz coordinates of each atom from xyz_file_name.xyz file
  [atom_count,atom_elem,atom_xyz] = get_xyz_data(xyz_file_name);
  %create color vector for atoms
  printf(" get color vector...");
  atom_cor = generate_color_vector (atom_count, atom_xyz, atom_elem);
  % read quaternions from xlsx_file_name.xlsx file [i j k r]
  printf(strcat(xlsx_file_name,"; "));
  Q(:,1:4) = xlsread (xlsx_file_name, task_name, "AN2:AQ9000" );
  % read xy gazepoint coordinates in screen fraction (needs to transform to pixels in config_pre_calc)
  printf("one more tab; ");
  gaze(:,1:2) = xlsread (xlsx_file_name, task_name, "D2:E9000" );
  printf("concluido!");toc();  %Time spent ~ 10 s
endif

if (config_pre_calc == 1)
  tic();printf("Calculando... ");
  ##encontrar centro de rotacao (boundingbox do jmol)
  atom_xyz(:,1:3) = normalize_jmol_rot_center (atom_xyz); %centralizado aqui!
  ## gerar matriz de rotacao pelo tempo com base nas coordenadas Q4 das rotacoes registradas
  for t = 1: size (Q,1) %criando matriz de rotacao para cada frame
    rot_vector(1:3,1:3,t) = rot_matrix (Q(t,4),Q(t,1),Q(t,2),Q(t,3));
    %rot_vector(1:3,1:3,1) = [0,1,0;-1,0,0;0,0,1];   %teste de rotacao em 90graus
    for a=1:atom_count     %para cada atomo, aplicar a matriz de rotacao
      %atom_xyzRot(a,1:3,t) = atom_xyz(a,1:3)*rot_vector(1:3,1:3,t); %1x3 * 3x3 = 1x3 (ESSE INVERTE A ROTACAO)
      atom_xyzRot(a,1:3,t) = (rot_vector(1:3,1:3,t)*atom_xyz(a,1:3)' )' ; %3x3 * 3x1 = 3x1 (ESSE E O CERTO) {rotacoes centralizadas no centro de rotacoes}
    endfor
  endfor
  for a=1:atom_count %gerar matriz de atomos do modelo referencia
    ref_atom_xyz(a,1:3) = (rot_matrix (config_ref_quat(1),config_ref_quat(2),config_ref_quat(3),config_ref_quat(4))*atom_xyz(a,1:3)')' ;
  endfor

  ##Gerar gaze coords. em pixel na tela (gaze_px), partindo do centro da referencia
  ##(gaze_ref_px) e do centro do canvas (gaze_cvs_px)
  gaze_px = gaze(:,1:2).*[config_screen_size];
  gaze_ref_px = gaze_px - config_ref_center_px;
  gaze_cvs_px(:,1:2) = gaze_px - config_cvs_center_px;

  ## gerar projecao xy dos atomos rotacionados pelas matrizes no tempo
  atom_xy(:,1:2,:) = atom_xyzRot(:,1:2,:); %atom_xy(atomos, xy, frame) {centralizado canvas}
  atom_xy_px(:,:,:) = config_fator_px*atom_xy(:,:,:); %obter pixels das coordenadas da projecao xy dos atomos {centralizado no centro de rot.}
  %ref_atom_xy(:,1:2) = ref_atom_xyz(:,1:2); %ref_atom_xy(atomo,xy) {centralizado na ref}
  ref_atom_xy_px(:,1:2) = ref_atom_xyz(:,1:2)*config_fator_px; %obter pixel xy
endif

## gerar vetor convex hull no tempo
if (config_convhull == 1)
  for t = 1:size(Q,1) %aplicar convex hull em cada frame.
    [H,area(t,1)] = convhull (atom_xy(:,1,t),atom_xy(:,2,t)); %se deixar [H,v]= ..., calcula volume em V
  endfor
  area(:,2) = atom_count./area(:,1);   %densidade de atomos na area do convexhull
  if ( isempty( strfind(config_xls_write,"1") ) != true)
    printf("desativei esta funcao xlswrite");
    %xlswrite (strcat("compiladoPorSubjects_",task(config_task),".xlsx"), area, num2str(subject), "AD2");   %registrar area do convexhull e densidade de �tomos nessa �rea
  endif
  printf("...convhull ok.");
endif

## gerar distancias do gazepoint(fracao da tela) pra cada atomo
if (config_gaze_atom_dist == 1)
  %gaze_atom_dist é = sqrt((x-x')^2 + (y-y')^2 )
  gaze_ref_atom_dist = gaze_atom_dist = zeros (size(Q,1),atom_count);
  for a=1:atom_count
    for t=1 : size(Q,1)
      gaze_atom_dist(t,a) = sqrt ( (atom_xy_px(a,1,t)-gaze_cvs_px(t,1))^2 + (atom_xy_px(a,2,t)-gaze_cvs_px(t,2))^2 );
      gaze_ref_atom_dist(t,a) = sqrt ( (ref_atom_xy_px(a,1)-gaze_cvs_px(t,1))^2 + (ref_atom_xy_px(a,2)-gaze_cvs_px(t,2))^2 );
    endfor
  endfor
  for t=1 : size(Q,1) %varredura do tempo para preencher aonde esta o gaze e o atomo mais proximo do gaze
    if ( -200<gaze_cvs_px(t,1) && gaze_cvs_px(t,1)<200 && -200<gaze_cvs_px(t,2) && gaze_cvs_px(t,2)<200  ) %marca em qual parte esta o gaze da pessoa
      gaze_status(t,1) = 2;
      [atom_closer_to_gaze(t,1),atom_closer_to_gaze(t,2)] = min (gaze_atom_dist(t,:));
    elseif ( -604<gaze_cvs_px(t,1) && gaze_cvs_px(t,1)<-204 && -203<gaze_cvs_px(t,2) && gaze_cvs_px(t,2)<197  )
      gaze_status(t,1) = 1;
      [atom_closer_to_gaze(t,1),atom_closer_to_gaze(t,2)] = min (gaze_ref_atom_dist(t,:));
    else
      gaze_status(t,1) = 0;
      atom_closer_to_gaze(t,1:2) = [0,0];
    endif
  endfor
  printf("...gaze dist ok.");
  if ( isempty( strfind(config_xls_write,"2") ) != true)
    xlswrite (xlsx_file_name, [gaze_cvs_px,gaze_status,atom_closer_to_gaze], num2str(subject), "AV2");   %registrar area do convexhull e densidade de atomos nessa area
  endif
  printf("...e impresso.");
printf("calculos concluidos! ");toc(); %Time spent (gaze dist)~ 3 s
endif

## calcular gradiende de transparencia da molecula baseado no tempo de proximidade do gaze
if (config_dist_integral == 1)
  printf("Cruzando dados do gaze com a posicao dos atomos...");tic();
  atom_gaze_alfa = ref_atom_gaze_alfa = zeros (atom_count,1);
  ref_dist = dist = zeros (size (Q,1),2,atom_count);
  ref_dist_exp = dist_exp = zeros (size(Q,1), atom_count);
  for a=1 : atom_count
    ##primeiro atomos da referencia : integral ( exp( - ( (x(t)-cx)^2 + (y(t)-cy)^2)/ (2*config_gauss_wdt^2)) dt)
    for t=1 : size (Q,1)
      ref_dist(t,1:2,a) = [ (gaze_ref_px(t,1) - ref_atom_xy_px(a,1)).^2 , (gaze_ref_px(t,2) - ref_atom_xy_px(a,2)).^2 ];
      dist(t,1:2,a) = [ (gaze_cvs_px(t,1) - atom_xy_px(a,1,t)).^2 , (gaze_cvs_px(t,2) - atom_xy_px(a,2,t)).^2 ];
    endfor
  endfor
  for c=1:10
    for a=1 : atom_count
      config_gauss_wdt = c*20;
      ref_dist_exp = exp (-(ref_dist(:,1,a)+ref_dist(:,2,a))/(2*config_gauss_wdt^2) ) ; %calcula e multiplica por 0/1 logico, se gaze esta dentro da ref
      dist_exp = exp (-(dist(:,1,a)+dist(:,2,a))/(2*config_gauss_wdt^2) ) ; %idem, mas pra referencia
      ref_atom_gaze_alfa(a) = sum(ref_dist_exp.*(gaze_status==1) ); %soma tudo
      atom_gaze_alfa(a) = sum(dist_exp.*(gaze_status==2) );
    endfor
    plot(atom_gaze_alfa);
    xlswrite ("lista_alfas_atomos.xlsx", [config_gauss_wdt;0;ref_atom_gaze_alfa], strcat("ref_",task_name), strcat (char (c+64), "2") );
    xlswrite ("lista_alfas_atomos.xlsx", [config_gauss_wdt;0;atom_gaze_alfa], strcat(task_name), strcat (char (c+64), "2") );
  endfor
  printf("concluido!");toc(); %Time spent (gaze dist)~ 60 s
endif

## calcular gradiende de transparencia da molecula VARIANDO NO TEMPO
if (config_dist_integral_temporal == 1)
  printf("Gerando animacao do gradiente de transparencia...");tic();
  atom_gaze_alfa = ref_atom_gaze_alfa = zeros (atom_count,1); #(a,1)
  ref_dist = dist = zeros (size (Q,1),2,atom_count);  #(t,1:2,a)
  ref_dist_exp = dist_exp = zeros (size (Q,1), atom_count); #(t,a)
  ##formula da gaussiana: integral ( exp( - ( (x(t)-cx)^2 + (y(t)-cy)^2)/ (2*config_gauss_wdt^2)) dt)
  for a=1 : atom_count
    for t=1 : size (Q,1)
      ref_dist(t,1:2,a) = [ (gaze_ref_px(t,1) - ref_atom_xy_px(a,1)).^2 , (gaze_ref_px(t,2) - ref_atom_xy_px(a,2)).^2 ];
      dist(t,1:2,a) = [ (gaze_cvs_px(t,1) - atom_xy_px(a,1,t)).^2 , (gaze_cvs_px(t,2) - atom_xy_px(a,2,t)).^2 ];
    endfor
  endfor
%  for c=1:10
  for c=1:1
    for a=1 : atom_count %preenchendo valores das integrais da gaussiana
      config_gauss_wdt = c*20; %largura da gaussiana
      ref_dist_exp(:,a) = exp (-(ref_dist(:,1,a)+ref_dist(:,2,a))/(2*config_gauss_wdt^2) ) ; %calcula e multiplica por 0/1 logico, se gaze esta dentro da ref
      dist_exp(:,a) = exp (-(dist(:,1,a)+dist(:,2,a))/(2*config_gauss_wdt^2) ) ; %idem, mas pra referencia
    endfor
    janela = round (size (Q,1)*0.1); %tamanho da janela de calculo do gradiente
    for t=1 : size (Q,1) %montando tabela do gradiente de transparencia com janela de frames
      s = max([t-janela, 1]); %regra pra janela(s:t)
      ref_atom_gaze_alfa(t,1:atom_count) = sum (ref_dist_exp(s:t,1:atom_count).*(gaze_status(s:t)==1) ); %soma tudo e armazena pra um atomo
      atom_gaze_alfa(t,1:atom_count) = sum (dist_exp(s:t,1:atom_count).*(gaze_status(s:t)==2) );
        %INSERIR NORMALIZACAO AQUI
    endfor
    plot(atom_gaze_alfa(end,:));
    xlswrite ("teste.xlsx", [ref_atom_gaze_alfa], strcat("ref_",task_name,"_",num2str (config_gauss_wdt),"_RAW"), "E2" );
    xlswrite ("teste.xlsx", [Q,atom_gaze_alfa], strcat(task_name,"_",num2str (config_gauss_wdt),"_RAW"), "A2" );
%    xlswrite ("lista_alfas_atomos.xlsx", [ref_atom_gaze_alfa], strcat("ref_",task_name,"_",num2str (config_gauss_wdt)), "E2" );
%    xlswrite ("lista_alfas_atomos.xlsx", [Q,atom_gaze_alfa], strcat(task_name,"_",num2str (config_gauss_wdt)), "A2" );
  endfor
  printf("concluido!");toc();
endif

if (2==1) %testes de sigma
  for a=1 : atom_count
    config_gauss_wdt = 10;
    ref_dist_exp = exp (-(ref_dist(:,1,a)+ref_dist(:,2,a))/(2*config_gauss_wdt^2) ) ; %calcula e multiplica por 0/1 logico, se gaze esta dentro da ref
    ref_atom_gaze_alfa(a) = sum(ref_dist_exp.*(gaze_status==1) );
    dist_exp = exp (-(dist(:,1,a)+dist(:,2,a))/(2*config_gauss_wdt^2) ) ; %calcula e multiplica por 0/1 logico, se gaze esta dentro da ref
    atom_gaze_alfa(a) = sum(dist_exp.*(gaze_status==1) );
  endfor
  xlswrite ("lista_alfas_atomos.xlsx", [config_gauss_wdt;0;ref_atom_gaze_alfa], strcat("ref_",task_name), strcat ("K", "2") );
  xlswrite ("lista_alfas_atomos.xlsx", [config_gauss_wdt;0;atom_gaze_alfa], task_name, strcat ("K", "2") );
endif

## graficos
if ( isempty( strfind(config_graficos,"1") ) != 1) %scatter 2D da molecula sem rotacao
  figure (1); %proj xy
    scatter (atom_xyz(:,1),atom_xyz(:,2), atom_cor(:,1), atom_cor(:,2:4));
    title ("Projecao xy dos vertices");
    axis ("square", "equal");
    xlabel ("x"); ylabel ("y");
endif
if ( isempty( strfind(config_graficos,"2") ) != 1) %scatter 2D da projecao xy ROTACIONADA no instante frame
  frame=1;
  figure (2);
%    plot (atom_xy(:,1,frame)(H), atom_xy(:,2,frame)(H), "r-");  %contorno
    scatter (atom_xy(:,1,frame),atom_xy(:,2,frame), atom_cor(:,1), atom_cor(:,2:4));
    title ("Projecao xy rotacionada dos vertices");
    axis ("square", "equal");
    xlabel ("x"); ylabel ("y");
endif
if ( isempty( strfind(config_graficos,"3") ) != 1) %scatter 3D da molecula sem rotacao
  figure (3);
    scatter3 (atom_xyz(:,1), atom_xyz(:,2), atom_xyz(:,3), atom_cor(:,1), atom_cor(:,2:4));
    title ("Scatter 3D dos vertices");
    axis ("square", "equal");
    xlabel("x"); ylabel("y"); zlabel("z");
endif
if ( isempty( strfind(config_graficos,"4") ) != 1) %scatter 3D da molecula rotacionada no instante frame
  figure (4);
    frame = 1;
    scatter3 (atom_xyzRot(:,1,frame), atom_xyzRot(:,2,frame), atom_xyzRot(:,3,frame), atom_cor(:,1), atom_cor(:,2:4));
    title ("Scatter 3D dos vertices com a orientação em t=frame");
    axis ("square", "equal");
    xlabel("x"); ylabel("y"); zlabel("z");
endif
if ( isempty( strfind(config_graficos,"5") ) != 1) %scatter 3D da molecula rotacionada no instante frame
  figure (5);
%    for (frame=1 : size(Q,1) )
    for (frame=100 : 110 )
      atom_cor_proximidade_gaze = atom_cor;
      if ( gaze_status(frame,1) == 2) %se gaze esta no canvas 2, destacar atomo mais proximo
        temp_atom_cor = atom_cor_proximidade_gaze(atom_closer_to_gaze(frame,2),:);
        temp_atom_cor = [temp_atom_cor(1,1)*3, 0,1,0]; %triplica o tamanho do atomo mais proximo do gaze e fica verde
        atom_cor_proximidade_gaze(atom_closer_to_gaze(frame,2),:) = temp_atom_cor;
      endif
      scatter (atom_xy_px(:,1,frame), atom_xy_px(:,2,frame), atom_cor_proximidade_gaze(:,1), atom_cor_proximidade_gaze(:,2:4));
      pause(0.2);
    endfor
    axis ([-200,200,-200,200]);
    xlabel("x"); ylabel("y");
endif

