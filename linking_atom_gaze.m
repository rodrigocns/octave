pkg load io
pkg load quaternion

%{
 v0.9 - get ready for the manuscript
 Dados da tarefa (task) agora estao dentro do arquivo de funcao get_task_data.m
%}
%{
NOTES: atoms and vertices are used interchangably here because each atom is
basically a vertice in 3D space of a known element contained in the .xyz file.
This "element" dictates its size and color inside the generate_color_vector()
function for any plot to be rendered at the end of the script, but it is not used in any calculations inside this script
%}

% ===CONFIGS===
task = ["C","G","I"];      %identifiers of each task. Rename accordingly.
config_task_index =     2; %Choose what task will be used from "task" array, above
config_input_data =     1; %1 to read .xlsx input files
config_pre_calc =       1; %1 to execute precalculations (needed for most functions)
config_convhull =       0; %1 to calculate convexhull of the atoms xy-projection in each time frame
config_gaze_atom_dist = 1; %1 to calculate distance matrix between gazepoint and interactive atoms/vertices.
config_dist_integral =  1; %1 to calculate transparency array of each atom (least seen, most transparent)
config_dist_integral_temporal =  1; %1 to calculate transparency array of each atom WITH TIME (least seen, most transparent)
config_xls_write =  "0"; %1: convexhull, 2:gazepoint. String that contains which files to write. (if string contains the char "1" it will write the convexhull data in a file)
config_graficos =   "0"; %1:2D Scatter. 2:2D Scatter in time, 3:3D Scatter,  4:3D Scatter in time, 5:(DEBUG) 3D scatter in time with nearest gaze
subject =  63; %1~62
config_screen_size = [1360,720];
config_ref_center_px = [211,376];
config_cvs_center_px = [615,379];

config_factor_px = 23.699; %screen-fraction x pixel-position ratio (formulas in User 1_all_gaze.xlsx > "px to x coord")
config_gauss_wdt = 0.922; %std deviation of config_factor_px

%xlsx_file_name = strcat ("compiladoPorSubjects_", task(config_task_index), ".xlsx");
xlsx_file_name = strcat ("tabela_de_dados_opengaze",".xlsx");
[config_model_name,task_name,config_ref_quat] = get_task_data (config_task_index);
xyz_file_name = strcat ("modelos/", config_model_name, ".xyz"); %name of .xyz file used in task

% ===FUNCTIONS===

% returns toration matrix {R} from a quaternion {qr,qi,qj,qk}
function R = rot_matrix (qr,qi,qj,qk, s = 1)
  %dont forget jmol data uses this order: qi,qj,qk,qr.
  R= [1-2*s*(qj^2 + qk^2), 2*s*(qi*qj - qk*qr), 2*s*(qi*qk + qj*qr);
      2*s*(qi*qj + qk*qr), 1-2*s*(qi^2 + qk^2), 2*s*(qj*qk - qi*qr);
      2*s*(qi*qk - qj*qr), 2*s*(qj*qk + qi*qr), 1-2*s*(qi^2 + qj^2)];
endfunction

% codify atoms by element|| in:(atom_count,atom_xyz,atom_elem) || out:[size,R,G,B]
function atom_cor = generate_color_vector (atom_count, atom_xyz, atom_elem)
  %used in some graph renderings at the end of this script
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

% normalize xyz coords so its center of rotation (center of jmol boundingbox ) is 0,0,0
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
  %take [number_id, element, xyz coordinates] of each atom from xyz_file_name.xyz file
  [atom_count,atom_elem,atom_xyz] = get_xyz_data(xyz_file_name);
  %create color vector for atoms
  printf(" get color vector...");
  atom_cor = generate_color_vector (atom_count, atom_xyz, atom_elem);
  % read quaternions from xlsx_file_name.xlsx file [i j k r]
  printf(strcat(xlsx_file_name,"; "));
  Q(:,1:4) = xlsread (xlsx_file_name, task_name, "AN2:AQ9000" );
  % read xy gazepoint coordinates in screen fraction (needs to transform to pixels in config_pre_calc)
  printf("one more tab... ");
  gaze(:,1:2) = xlsread (xlsx_file_name, task_name, "D2:E9000" );
  printf("done! ");toc();  %Time spent ~ 10 s
endif

% pre-calculation of things
if (config_pre_calc == 1)
  tic();printf("Calculating... ");
  % find rotation center (jmol boundingbox)
  atom_xyz(:,1:3) = normalize_jmol_rot_center (atom_xyz); %centralized here!
  % create array of rotation matrix in time from quaternions
  for t = 1: size (Q,1) %create rotation matrix for each frame
    %rot_vector(1:3,1:3,1) = [0,1,0;-1,0,0;0,0,1];   %DEBUG: this line should rotate in 90degrees
    rot_vector(1:3,1:3,t) = rot_matrix (Q(t,4),Q(t,1),Q(t,2),Q(t,3));
    for a=1:atom_count     %apply rotation for each atom.
      %3x3 * 3x1 = 3x1 {rotation center at 0,0,0}
      % this is the right order. changing it will give reversed results!
      atom_xyzRot(a,1:3,t) = (rot_vector(1:3,1:3,t)*atom_xyz(a,1:3)' )' ;
    endfor
  endfor
  %create atom matrix of the reference model
  for a=1:atom_count
    ref_atom_xyz(a,1:3) = (rot_matrix (config_ref_quat(1),config_ref_quat(2),config_ref_quat(3),config_ref_quat(4))*atom_xyz(a,1:3)')' ;
  endfor

  % Generate gaze coordinates in screen px (gaze_px), from reference center
  %(gaze_ref_px) and from canvas center (gaze_cvs_px)
  gaze_px = gaze(:,1:2).*[config_screen_size];
  gaze_ref_px = gaze_px - config_ref_center_px;
  gaze_cvs_px(:,1:2) = gaze_px - config_cvs_center_px;

  % Create temporal xy projection of atoms rotated by matrices
  atom_xy(:,1:2,:) = atom_xyzRot(:,1:2,:); %atom_xy(atoms, xy, frame) {centralized canvas}
  atom_xy_px(:,:,:) = config_factor_px*atom_xy(:,:,:); %get px coordinates of atoms xy projection {already centralized}
  %ref_atom_xy(:,1:2) = ref_atom_xyz(:,1:2); %ref_atom_xy(atom,xy) {centralized in ref}
  ref_atom_xy_px(:,1:2) = ref_atom_xyz(:,1:2)*config_factor_px; %get pixel xy
endif

% Create convex hull time array
if (config_convhull == 1)
  %apply convex hull in each frame.
  for t = 1:size(Q,1)
    [H,area(t,1)] = convhull (atom_xy(:,1,t),atom_xy(:,2,t));
  endfor
  %atom density in convexhull area
  area(:,2) = atom_count./area(:,1);
  if ( isempty( strfind(config_xls_write,"1") ) != true)
    %write file with convexhull area and atom density in said area
    %xlswrite (strcat("compiladoPorSubjects_",task(config_task_index),".xlsx"), area, num2str(subject), "AD2");
    printf("DEBUG: xlswrite was turned off");
  endif
  printf("...convhull ok.");
endif

% Generate distances from gazepoint(screen fraction) for each atom
if (config_gaze_atom_dist == 1)
  %initialize arrays (time x atom_count)
  printf("Calculating gaze-atom distance array.."); tic();
  gaze_ref_atom_dist = gaze_atom_dist = zeros (size(Q,1),atom_count);
  for a=1:atom_count %for each atom
    for t=1 : size(Q,1) %for each point in time
      %Formula: gaze_atom_dist = sqrt( (x-x')^2 + (y-y')^2 )
      gaze_atom_dist(t,a) = sqrt ( (atom_xy_px(a,1,t)-gaze_cvs_px(t,1))^2 + (atom_xy_px(a,2,t)-gaze_cvs_px(t,2))^2 );
      gaze_ref_atom_dist(t,a) = sqrt ( (ref_atom_xy_px(a,1)-gaze_cvs_px(t,1))^2 + (ref_atom_xy_px(a,2)-gaze_cvs_px(t,2))^2 );
    endfor
  endfor
  %time scan to fill in where are the gaze and its nearest atom
  for t=1 : size(Q,1)
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
  printf(".array calculated.");
  if ( isempty( strfind(config_xls_write,"2") ) != true)
    %write convexhull area and atom density from these columns
    xlswrite (xlsx_file_name, [gaze_cvs_px,gaze_status,atom_closer_to_gaze], num2str(subject), "AV2");
    printf(".file written.");
  endif
  printf(".complete! ");toc(); %Time spent (gaze dist)~ 3 s
endif

% Calculate 3D object transparency gradient from time spent in proximity of gaze
if (config_dist_integral == 1)
  printf("Calculating transparency gradient.");tic();
  %initialize
  atom_gaze_alfa = ref_atom_gaze_alfa = zeros (atom_count,1);
  ref_dist = dist = zeros (size (Q,1),2,atom_count);
  ref_dist_exp = dist_exp = zeros (size(Q,1), atom_count);
  for a=1 : atom_count
    %gaussian formula: integral ( exp( - ( (x(t)-cx)^2 + (y(t)-cy)^2)/ (2*config_gauss_wdt^2)) dt)
    for t=1 : size (Q,1)
      ref_dist(t,1:2,a) = [ (gaze_ref_px(t,1) - ref_atom_xy_px(a,1)).^2 , (gaze_ref_px(t,2) - ref_atom_xy_px(a,2)).^2 ];
      dist(t,1:2,a) = [ (gaze_cvs_px(t,1) - atom_xy_px(a,1,t)).^2 , (gaze_cvs_px(t,2) - atom_xy_px(a,2,t)).^2 ];
    endfor
  endfor
  for c=1:10
    for a=1 : atom_count
      config_gauss_wdt = c*20;
      %calcula e multiplica por 0/1 logico, se gaze esta dentro da ref
      % calculate and multiply by binary (0/1), testing for gaze is inside reference
      ref_dist_exp = exp (-(ref_dist(:,1,a)+ref_dist(:,2,a))/(2*config_gauss_wdt^2) ) ;

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

% Calculate 3D object transparency gradient from time spent in proximity of gaze
if (config_dist_integral_temporal == 1)
  printf("Calculating transparency gradient animation.");tic();
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

if (2==1) %sigma tests
  for a=1 : atom_count
    config_gauss_wdt = 10;
    %calcula e multiplica por 0/1 logico, se gaze esta dentro da ref
    ref_dist_exp = exp (-(ref_dist(:,1,a)+ref_dist(:,2,a))/(2*config_gauss_wdt^2) ) ;
    ref_atom_gaze_alfa(a) = sum(ref_dist_exp.*(gaze_status==1) );
    %calcula e multiplica por 0/1 logico, se gaze esta dentro da ref
    dist_exp = exp (-(dist(:,1,a)+dist(:,2,a))/(2*config_gauss_wdt^2) ) ;
    atom_gaze_alfa(a) = sum(dist_exp.*(gaze_status==1) );
  endfor
  xlswrite ("lista_alfas_atomos.xlsx", [config_gauss_wdt;0;ref_atom_gaze_alfa], strcat("ref_",task_name), strcat ("K", "2") );
  xlswrite ("lista_alfas_atomos.xlsx", [config_gauss_wdt;0;atom_gaze_alfa], task_name, strcat ("K", "2") );
endif

% plots

if ( isempty( strfind(config_graficos,"1") ) != 1) %2D scatter of vertices/atoms without rotation
  figure (1); %proj xy
    scatter (atom_xyz(:,1),atom_xyz(:,2), atom_cor(:,1), atom_cor(:,2:4));
    title ("xy-projection of vertices");
    axis ("square", "equal");
    xlabel ("x"); ylabel ("y");
endif
if ( isempty( strfind(config_graficos,"2") ) != 1) %2D scatter of vertices/atoms rotated in t= "frame"
  frame=1;
  figure (2);
%    plot (atom_xy(:,1,frame)(H), atom_xy(:,2,frame)(H), "r-");  %border
    scatter (atom_xy(:,1,frame),atom_xy(:,2,frame), atom_cor(:,1), atom_cor(:,2:4));
    title (strcat("Rotated xy-projection of vertices in frame ", num2str(frame) ) );
    axis ("square", "equal");
    xlabel ("x"); ylabel ("y");
endif
if ( isempty( strfind(config_graficos,"3") ) != 1) %3D scatter of vertices/atoms without rotation
  figure (3);
    scatter3 (atom_xyz(:,1), atom_xyz(:,2), atom_xyz(:,3), atom_cor(:,1), atom_cor(:,2:4));
    title ("3D Scatter of vertices");
    axis ("square", "equal");
    xlabel("x"); ylabel("y"); zlabel("z");
endif
if ( isempty( strfind(config_graficos,"4") ) != 1) %3D scatter of vertices/atoms rotated in t= "frame"
  figure (4);
    frame = 1;
    scatter3 (atom_xyzRot(:,1,frame), atom_xyzRot(:,2,frame), atom_xyzRot(:,3,frame), atom_cor(:,1), atom_cor(:,2:4));
    title (strcat("Rotated 3D Scatter of vertices in frame ", num2str(frame) ) );
    axis ("square", "equal");
    xlabel("x"); ylabel("y"); zlabel("z");
endif
if ( isempty( strfind(config_graficos,"5") ) != 1) %3D scatter of rotated object from time index "f_0" to "f_end"
  %TBD
  figure (5);
    f_0 = 100;
    f_end = 110;
    %use 1 : size(Q,1) for the entire task
    for (frame=f_0 : f_end )
      atom_color_gaze_proximity = atom_cor;
      %only highlight closest atom/point if gaze is inside the interactive window
      % to avoid visual artifacts related to the reference obj. window located
      % at the left side
      if ( gaze_status(frame,1) == 2)
        temp_atom_cor = atom_color_gaze_proximity(atom_closer_to_gaze(frame,2),:);
        %triples atom size closest to gaze and paint it green
        temp_atom_cor = [temp_atom_cor(1,1)*3, 0,1,0];
        atom_color_gaze_proximity(atom_closer_to_gaze(frame,2),:) = temp_atom_cor;
      endif
      scatter (atom_xy_px(:,1,frame), atom_xy_px(:,2,frame), atom_color_gaze_proximity(:,1), atom_color_gaze_proximity(:,2:4));
      pause(0.2);
    endfor
    axis ([-200,200,-200,200]);
    xlabel("x"); ylabel("y");
endif

