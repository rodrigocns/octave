pkg load io
pkg load quaternion
#http://jmol.sourceforge.net/demo/jssample0.html

function R = rot_matrix (qr,qi,qj,qk, s = 1) #no .xlsx está qi,qj,qk,qr
  R= [1-2*s*(qj^2 + qk^2), 2*s*(qi*qj - qk*qr), 2*s*(qi*qk + qj*qr);
      2*s*(qi*qj + qk*qr), 1-2*s*(qi^2 + qk^2), 2*s*(qj*qk - qi*qr);
      2*s*(qi*qk - qj*qr), 2*s*(qj*qk + qi*qr), 1-2*s*(qi^2 + qj^2)]; 
endfunction

function Q = axangle2quat (x,y,z,angle) #Q = [qw,qx,qy,qz]. angle in degrees
  Q = [ cos(deg2rad(angle/2)), x*sin(deg2rad(angle/2)), y*sin(deg2rad(angle/2)), z*sin(deg2rad(angle/2))];
endfunction

function atom_cor = gerar_vetor_cores (atom_count, atom_xyz, atom_elem)  #codificar átomos por elemento|| in:(atom_count,atom_xyz,atom_elem) || out:[size,R,G,B]
  for i = 1:atom_count
    switch ( strvcat(atom_elem(i)) )  #strvcat extrai a string de um cell array
      case "H"
        atom_xyz(i,4) = 1;
        atom_cor(i,1:4) = [74, 0.75,0.75,0.75];
      case "C"
        atom_xyz(i,4) = 12;
        atom_cor(i,1:4) = [154, 0,0,0];
      case "N"
        atom_xyz(i,4) = 14;
        atom_cor(i,1:4) = [140, 0,0,1];
      case "O"
        atom_xyz(i,4) = 16;
        atom_cor(i,1:4) = [146, 1,0,0];
      otherwise
        error ("ERRO! Novo elemento detectado. Alterar código");
        atom_xyz(i,4) = -1;
        atom_cor(i,1:4) = [20, 0,1,0];
    endswitch
  endfor
endfunction

function norm_atom_xyz = normalize_jmol_rot_center (atom_xyz) #normaliza coordenadas de amtomos da matriz para o centro de rotacao ser 0,0,0
  max_xyz = max(atom_xyz(:,1:3));
  min_xyz = min(atom_xyz(:,1:3));
  correcao_centro = (max_xyz+min_xyz)/2;
  norm_atom_xyz = atom_xyz(:,1:3) - correcao_centro;
endfunction


## ===CONFIGS===
tipo = ["C","G","I"];
config_tipo =           2; #1:C, 2:G, 3:I
config_input_data =     0; #1 se vai limpar dados e ler dos .xlsx
config_pre_calc =       0; #1 se for realizar outros calculos (importantes para o resto)
config_convhull =       0; #1 para calcular a convexhull da projeção xy dos atomos em cada tempo
config_gaze_atom_dist = 0; #1 para calcular matriz de distancia entre gazepoint e atomos interativos.
config_dist_integral =  1; #1 para calcular coluna de alfa/transparencia de cada atomo menos visto
config_xls_write =  "0"; #1: convexhull, 2:gazepoint
config_graficos =   "0"; #1:projecao, 2:projecao em t, 3:scatter3,  4:scatter3 em frame, 5:scatter3 em frame com gaze mais próximo
subject =  63; #1~62
config_screen_size = [1360,720];
config_ref_center_px = [211,376];
config_cvs_center_px = [615,379];
config_fator_px = 23.699; #razao entre fracao da tela e posicao em px (calc. no xlsx)
config_gauss_wdt = 0.922; #desvio padrao do config_fator_px
#xlsx_file_name = strcat ("compiladoPorSubjects_", tipo(config_tipo), ".xlsx");
xlsx_file_name = strcat ("tabela_de_dados_opengaze.xlsx");
xyz_file_name = "convexhull/xyz coords.xlsx";
## ===PREENCHER===
if (config_tipo == 1 ) #nome do modelo molecular
  model_name = 'pseudobatracotoxin_molecule'; #caprolactama, pseudobatracotoxin_molecule
  tipo_name = 'cor';
  config_ref_q = [0.6873,  -0.4903,  -0.2825,  -0.4554]; #quaternios para a referencia
elseif (config_tipo == 2)
  model_name = 'pseudobatracotoxin_molecule';
  tipo_name = 'cinza';
  config_ref_q = [0.6873,  -0.4903,  -0.2825,  -0.4554]; #quaternios para a referencia
elseif (config_tipo == 3)
  model_name = 'caprolactama';
  tipo_name = 'intro';
  config_ref_q = [0.5183,   0.5105,   0.4575,  -0.5105]; #quaternios para a referencia
else
  error ("ERRO! Tipo de teste não reconhecido!");
endif
## ==============

if (config_input_data == 1)
  clear -x config* subject *name 
  tic();printf("Lendo arquivos xlsx... ");
  ##pegar contagem, lista de coordenadas xyz e tipo de cada átomo do modelo.
  atom_count = xlsread (xyz_file_name, model_name, 'A1:A1');
  printf("1.");
  [atom_xyz, atom_elem] = xlsread (xyz_file_name, model_name, strcat('A3:D',num2str(atom_count+2) ));
  printf("2.");
  atom_cor = gerar_vetor_cores (atom_count, atom_xyz, atom_elem); #gerar vetor de cores para os atomos
  ## ler quaternios do .xlsx [x y z theta]
  Q(:,1:4) = xlsread (xlsx_file_name, tipo_name, "AN2:AQ9000" );
  printf("3.");
  ## ler coordenadas x,y do gazepoint (está em fracao da tela, precisa corrigir para pixels)
  gaze(:,1:2) = xlsread (xlsx_file_name, tipo_name, "D2:E9000" );
  printf("4.");
  printf("concluído!");toc();
endif

if (config_pre_calc == 1)
  tic();printf("Calculando... ");
  ##encontrar centro de rotação (boundingbox do jmol)
  atom_xyz(:,1:3) = normalize_jmol_rot_center (atom_xyz); #centralizado aqui!
  ## gerar matriz de rotação pelo tempo com base nas coordenadas Q4 das rotações registradas
  for t = 1: size (Q,1) #para cada frame\tempo
    rot_vector(1:3,1:3,t) = rot_matrix (Q(t,4),Q(t,1),Q(t,2),Q(t,3));  #criando a matriz de rotação
    #rot_vector(1:3,1:3,1) = [0,1,0;-1,0,0;0,0,1];   #teste de rotação em 90graus
    for a=1:atom_count     #para cada atomo, aplicar a matriz de rotação
      #atom_xyzRot(a,1:3,t) = atom_xyz(a,1:3)*rot_vector(1:3,1:3,t); #1x3 * 3x3 = 1x3 (ESSE INVERTE A ROTACAO)
      atom_xyzRot(a,1:3,t) = (rot_vector(1:3,1:3,t)*atom_xyz(a,1:3)' )' ; #3x3 * 3x1 = 3x1 ESSE E O CERTO) {rotacoes centralizadas no centro de rotacoes}
    endfor
  endfor
  for a=1:atom_count #gerar matriz de atomos do modelo referencia
    ref_atom_xyz(a,1:3) = (rot_matrix (config_ref_q(1),config_ref_q(2),config_ref_q(3),config_ref_q(4))*atom_xyz(a,1:3)')' ; 
  endfor

  ##gaze coords. em pixel na tela (gaze_px), partindo do centro da referencia 
  ##(gaze_ref_px) e do centro do canvas (gaze_cvs_px)
  gaze_px = gaze(:,1:2).*[config_screen_size];
  gaze_ref_px = gaze_px - config_ref_center_px;
  gaze_cvs_px(:,1:2) = gaze_px - config_cvs_center_px; 

  
  ## gerar projeção xy dos átomos rotacionados pelas matrizes no tempo
  atom_xy(:,1:2,:) = atom_xyzRot(:,1:2,:); #atom_xy(atomos, xy, frame) {centralizado canvas}
  atom_xy_px(:,:,:) = config_fator_px*atom_xy(:,:,:); #obter pixels das coordenadas da projecao xy dos atomos {centralizado no centro de rot.}
  ref_atom_xy(:,1:2) = ref_atom_xyz(:,1:2); #ref_atom_xy(atomo,xy) {centralizado na ref}
  ref_atom_xy_px(:,1:2) = ref_atom_xy(:,1:2)*config_fator_px; #obter pixel xy
endif
  
## gerar vetor convex hull no tempo
if (config_convhull == 1)
  for t = 1:size(Q,1) #aplicar convex hull em cada frame.
    [H,area(t,1)] = convhull (atom_xy(:,1,t),atom_xy(:,2,t)); #se deixar [H,v]= ..., calcula volume em V
  endfor
  area(:,2) = atom_count./area(:,1);   #densidade de átomos na área do convexhull
  if ( isempty( strfind(config_xls_write,"1") ) != true)
    printf("desativei esta funcao xlswrite");
    #xlswrite (strcat("compiladoPorSubjects_",tipo(config_tipo),".xlsx"), area, num2str(subject), "AD2");   #registrar area do convexhull e densidade de átomos nessa área
  endif
  printf("...convhull ok.");
endif

## gerar distâncias do gazepoint(fracao da tela) pra cada átomo
if (config_gaze_atom_dist == 1)
  #gaze_atom_dist = sqrt((x-x')² + (y-y')² )
  gaze_ref_atom_dist = gaze_atom_dist = zeros (size(Q,1),atom_count);
  for a=1:atom_count
    for t=1 : size(Q,1)
      gaze_atom_dist(t,a) = sqrt ( (atom_xy_px(a,1,t)-gaze_cvs_px(t,1))^2 + (atom_xy_px(a,2,t)-gaze_cvs_px(t,2))^2 );
      gaze_ref_atom_dist(t,a) = sqrt ( (ref_atom_xy_px(a,1)-gaze_cvs_px(t,1))^2 + (ref_atom_xy_px(a,2)-gaze_cvs_px(t,2))^2 );
    endfor
  endfor
  for t=1 : size(Q,1) #varredura do tempo para preencher aonda está o gaze e o átomo mais próximo do gaze
    if ( -200<gaze_cvs_px(t,1) && gaze_cvs_px(t,1)<200 && -200<gaze_cvs_px(t,2) && gaze_cvs_px(t,2)<200  ) #marca em qual parte está o gaze da pessoa
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
    xlswrite (xlsx_file_name, [gaze_cvs_px,gaze_status,atom_closer_to_gaze], num2str(subject), "AV2");   #registrar area do convexhull e densidade de átomos nessa área
  endif
  printf("...e impresso.");
printf("calculos concluídos! ");toc();
endif

## calcular alfa de cada atomo baseado no tempo de proximidade do gaze
if (config_dist_integral == 1)
  atom_gaze_alfa = zeros (atom_count,1);
  dist = zeros (size (Q,1),2,atom_count);
  dist_exp = zeros (size(Q,1), atom_count);
  for a=1 : atom_count
    ##primeiro atomos da referencia : integral ( exp( - ( (x(t)-cx)^2 + (y(t)-cy)^2)/ (2*config_gauss_wdt^2)) dt)
    for t=1 : size (Q,1)
        dist(t,1:2,a) = [ (gaze_ref_px(t,1) - ref_atom_xy_px(a,1)).^2 , (gaze_ref_px(t,2) - ref_atom_xy_px(a,2)).^2 ];
#      if (gaze_status(t)==1) #se gaze esta na ref
#        dist_exp = exp (-(dist(:,1,a)+dist(:,2,a))/(2*config_gauss_wdt^2) );
#       dist_exp = exp (-(dist(:,1,a)+dist(:,2,a)) );
#      endif
    endfor
  endfor
  for c=1:10
    for a=1 : atom_count
      config_gauss_wdt = c*20;
      dist_exp = exp (-(dist(:,1,a)+dist(:,2,a))/(2*config_gauss_wdt^2) ) ; #calcula e multiplica por 0/1 logico, se gaze esta dentro da ref
      atom_gaze_alfa(a) = sum(dist_exp.*(gaze_status==1) );
    endfor
    plot(atom_gaze_alfa);
    xlswrite ("lista_alfas_atomos.xlsx", [config_gauss_wdt;0;atom_gaze_alfa], tipo_name, strcat (char (c+64), "2") );
  endfor
endif

if (2==1)
  for a=1 : atom_count
    config_gauss_wdt = 10;
    dist_exp = exp (-(dist(:,1,a)+dist(:,2,a))/(2*config_gauss_wdt^2) ) ; #calcula e multiplica por 0/1 logico, se gaze esta dentro da ref
    atom_gaze_alfa(a) = sum(dist_exp.*(gaze_status==1) );
  endfor
  xlswrite ("lista_alfas_atomos2.xlsx", [config_gauss_wdt;0;atom_gaze_alfa], tipo_name, strcat ("K", "2") );
endif

## graficos
if ( isempty( strfind(config_graficos,"1") ) != 1) #scatter 2D da molecula sem rotacao
  figure (1); #proj xy 
    scatter (atom_xyz(:,1),atom_xyz(:,2), atom_cor(:,1), atom_cor(:,2:4));
    title ("Projeção xy");
    axis ("square", "equal");
    xlabel ("x"); ylabel ("y");
endif
if ( isempty( strfind(config_graficos,"2") ) != 1) #scatter 2D da projecao xy ROTACIONADA no instante frame
  frame=1;
  figure (2); 
#    plot (atom_xy(:,1,frame)(H), atom_xy(:,2,frame)(H), "r-");  #contorno
    scatter (atom_xy(:,1,frame),atom_xy(:,2,frame), atom_cor(:,1), atom_cor(:,2:4));
    title ("Projeção xy rotacionada");
    axis ("square", "equal");
    xlabel ("x"); ylabel ("y");
endif
if ( isempty( strfind(config_graficos,"3") ) != 1) #scatter 3D da molecula sem rotacao
  figure (3); 
    scatter3 (atom_xyz(:,1), atom_xyz(:,2), atom_xyz(:,3), atom_cor(:,1), atom_cor(:,2:4));
    axis ("square", "equal");
    xlabel("x"); ylabel("y"); zlabel("z");
endif
if ( isempty( strfind(config_graficos,"4") ) != 1) #scatter 3D da molecula rotacionada no instante frame 
  figure (4); 
    frame = 1;
    scatter3 (atom_xyzRot(:,1,frame), atom_xyzRot(:,2,frame), atom_xyzRot(:,3,frame), atom_cor(:,1), atom_cor(:,2:4));
    axis ("square", "equal");
    xlabel("x"); ylabel("y"); zlabel("z");
endif
if ( isempty( strfind(config_graficos,"5") ) != 1) #scatter 3D da molecula rotacionada no instante frame 
  figure (5); 
#    for (frame=1 : size(Q,1) )
    for (frame=100 : 110 )
      atom_cor_proximidade_gaze = atom_cor;
      if ( gaze_status(frame,1) == 2) #se gaze esta no canvas 2, destacar atomo mais proximo
        temp_atom_cor = atom_cor_proximidade_gaze(atom_closer_to_gaze(frame,2),:);
        temp_atom_cor = [temp_atom_cor(1,1)*3, 0,1,0]; #triplica o tamanho do átomo mais próximo do gaze e fica verde
        atom_cor_proximidade_gaze(atom_closer_to_gaze(frame,2),:) = temp_atom_cor;
      endif
      scatter (atom_xy_px(:,1,frame), atom_xy_px(:,2,frame), atom_cor_proximidade_gaze(:,1), atom_cor_proximidade_gaze(:,2:4));
      pause(0.2);
    endfor
    axis ([-200,200,-200,200]);
    xlabel("x"); ylabel("y");
endif

