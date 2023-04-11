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

clear

## ===CONFIGS===
tipo = ["C","G","I"];
config_tipo =     1; #1:C, 2:G, 3:I
config_convhull =   0; #1 para calcular a convexhull da projeção xy dos atomos em cada tempo
config_gaze_atom_dist = 1; #1 para calcular matriz de distancia entre gazepoint e atomos interativos.
config_xls_write = "2"; #1: convexhull, 2:gazepoint
config_graficos = "5"; #1:projecao simples, 2:plot 3D, 3:plot 3D rotacionado em cada frame
subject =  37; #1~62
screen_size = [1360,720];
#xlsx_file_name = strcat ("compiladoPorSubjects_", tipo(config_tipo), ".xlsx");
xlsx_file_name = strcat ("tabela_de_dados_opengaze.xlsx");
xyz_file_name = "convexhull/xyz coords.xlsx";
## =============
## ===PREENCHER===
if (config_tipo == 1 ) #nome do modelo molecular
  model_name = 'pseudobatracotoxin_molecule'; #caprolactama, pseudobatracotoxin_molecule
  tipo_name = 'cor';
elseif (config_tipo == 2)
  model_name = 'pseudobatracotoxin_molecule';
  tipo_name = 'cinza';
elseif (config_tipo == 3)
  model_name = 'caprolactama';
else
  error ("ERRO! Tipo de teste não reconhecido!");
endif
## ==============

clear Q rot_vector atom_xyzRot xy_proj H area
##pegar contagem, lista de coordenadas xyz e tipo de cada átomo do modelo.
atom_count = xlsread (xyz_file_name, model_name, 'A1:A1');
[atom_xyz, atom_elem] = xlsread (xyz_file_name, model_name, strcat('A3:D',num2str(atom_count+2) ));
atom_cor = gerar_vetor_cores (atom_count, atom_xyz, atom_elem); #gerar vetor de cores para os atomos
## ler quaternios do .xlsx [x y z theta]
Q(:,1:4) = xlsread (xlsx_file_name, tipo_name, "AN2:AQ9000" );
## ler coordenadas x,y do gazepoint (está em fracao da tela, precisa corrigir para pixels)
gaze(:,1:2) = xlsread (xlsx_file_name, tipo_name, "D2:E9000" );

##encontrar centro de rotação (boundingbox do jmol)
atom_xyz(:,1:3) = normalize_jmol_rot_center (atom_xyz); #centralizado aqui!

## gerar matriz de rotação pelo tempo com base nas coordenadas Q4 das rotações registradas
for t = 1:size(Q,1) #para cada frame\tempo
  rot_vector(1:3,1:3,t) = rot_matrix (Q(t,4),Q(t,1),Q(t,2),Q(t,3));  #criando a matriz de rotação
  #rot_vector(1:3,1:3,1) = [0,1,0;-1,0,0;0,0,1];   #teste de rotação em 90graus
  for a=1:atom_count     #para cada atomo, aplicar a matriz de rotação
    #atom_xyzRot(a,1:3,t) = atom_xyz(a,1:3)*rot_vector(1:3,1:3,t); #1x3 * 3x3 = 1x3 (ESSE INVERTE A ROTACAO)
    atom_xyzRot(a,1:3,t) = (rot_vector(1:3,1:3,t)*atom_xyz(a,1:3)' )' ; #3x3 * 3x1 = 3x1 ESSE E O CERTO) {rotacoes centralizadas no centro de rotacoes}
  endfor
endfor

## gerar projeção xy dos átomos rotacionados pelas matrizes no tempo
xy_proj(:,1:2,:) = atom_xyzRot (:,1:2,:); #xy_proj(atomos, x e y, frame) {centralizado no centro de rotacao}

## gerar vetor convex hull no tempo
if (config_convhull == 1)
  for t = 1:size(Q,1) #aplicar convex hull em cada frame.
    [H,area(t,1)] = convhull (xy_proj(:,1,t),xy_proj(:,2,t)); #se deixar [H,v]= ..., calcula volume em V
  endfor
  area(:,2) = atom_count./area(:,1);   #densidade de átomos na área do convexhull
  if ( isempty( strfind(config_xls_write,"1") ) != 1)
    print("desativei esta funcao xlswrite");
    #xlswrite (strcat("compiladoPorSubjects_",tipo(config_tipo),".xlsx"), area, num2str(subject), "AD2");   #registrar area do convexhull e densidade de átomos nessa área
  endif
endif

## gerar distâncias do gazepoint(fracao da tela) pra cada átomo
if (config_gaze_atom_dist == 1)
  gaze_cnvs2_coords(:,1:2) = (gaze(:,1:2).*[screen_size])-[615,379]; #coord. xy na tela do gazepoint transformado para coordenadas xy centralizadas no centro do canvas interativo
  xy_proj_px(:,:,:) = 23.753*xy_proj(:,:,:); #obter pixels das coordenadas da projecao xy dos atomos {centralizado no centro de rot.}
  #gaze_atom_dist = sqrt((x-x')² + (y-y')² )
  for a=1:atom_count
    for t=1 : size(Q,1)
      gaze_atom_dist(t,a) = sqrt ( (xy_proj_px(a,1,t)-gaze_cnvs2_coords(t,1))^2 + (xy_proj_px(a,2,t)-gaze_cnvs2_coords(t,2))^2 );
    endfor
  endfor
  for t=1 : size(Q,1) #varredura do tempo para preencher aonda está o gaze e o átomo mais próximo do gaze
    if ( -200<gaze_cnvs2_coords(t,1) && gaze_cnvs2_coords(t,1)<200 && -200<gaze_cnvs2_coords(t,2) && gaze_cnvs2_coords(t,2)<200  ) #marca em qual parte está o gaze da pessoa
      gaze_status(t,1) = 2;
      [atom_closer_to_gaze(t,1),atom_closer_to_gaze(t,2)] = min (gaze_atom_dist(t,:));
    if ( -604<gaze_cnvs2_coords(t,1) && gaze_cnvs2_coords(t,1)<-204 && -203<gaze_cnvs2_coords(t,2) && gaze_cnvs2_coords(t,2)<197  )
      gaze_status(t,1) = 1;
      [atom_closer_to_gaze(t,1),atom_closer_to_gaze(t,2)] = 0;
    else
      gaze_status(t,1) = 0;
      [atom_closer_to_gaze(t,1),atom_closer_to_gaze(t,2)] = 0;
    endif
  endfor
  
  if ( isempty( strfind(config_xls_write,"2") ) != 2)
    xlswrite (xlsx_file_name, [gaze_cnvs2_coords,gaze_status,atom_closer_to_gaze], num2str(subject), "AV2");   #registrar area do convexhull e densidade de átomos nessa área
  endif
endif

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
#    plot (xy_proj(:,1,frame)(H), xy_proj(:,2,frame)(H), "r-");  #contorno
    scatter (xy_proj(:,1,frame),xy_proj(:,2,frame), atom_cor(:,1), atom_cor(:,2:4));
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
    axis ([-200,200,-200,200]);
    xlabel("x"); ylabel("y");
#    for (frame=1 : size(Q,1) )
    for (frame=1 : 20 )
      atom_cor_proximidade_gaze = atom_cor;
      if ( gaze_status(frame,1) == 2) #se gaze esta no canvas 2, destacar atomo mais proximo
        temp_atom_cor = atom_cor_proximidade_gaze(atom_closer_to_gaze(frame,2),:);
        temp_atom_cor = [temp_atom_cor(1,1)*3, 0,1,0]; #triplica o tamanho do átomo mais próximo do gaze e fica verde
        atom_cor_proximidade_gaze(atom_closer_to_gaze(frame,2),:) = temp_atom_cor;
      endif
      scatter (xy_proj_px(:,1,frame), xy_proj_px(:,2,frame), xy_proj_px(:,3,frame), atom_cor_proximidade_gaze(:,1), atom_cor_proximidade_gaze(:,2:4));
      pause(0.1);
    endfor
endif
