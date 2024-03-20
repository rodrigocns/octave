% Find data for computing the "quaternionic worm-path"
%from the "filename" excel spreadsheet file of unpackaged iRT data from gsheets

pkg load io
clear;


% To be chosen: 'session ID' and 'name of task'
%session_ID = 1682699553789; %volunteer1
session_ID = "1682707472090"; %volunteer2
task_ID = "bolaBastao_c"; % molecule with balls and sticks
%task_ID = "poligonFill"; % spatial polygons
%task_ID = "mrt" % third task of the iRT
filename = "iRT data.xlsx"; % default unpackaged data filename


tic();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUT

%session ID
prpt_sessionID = inputdlg ("Insert session ID value. Ex: 1682699553789, 1682707472090", "Input Session ID");
if (isempty ( char (prpt_sessionID) ) == 1 )
  printf ( cstrcat ("Blank input, using default value: ",session_ID,"\n"));
else
  session_ID = char(prpt_sessionID);
endif

%task ID
prpt_taskID = inputdlg ("Insert task ID value. Ex: bolaBastao\_c, poligonFill, mrt", "Input Task ID");
if (isempty ( char (prpt_taskID) ) == 1 )
  printf ( cstrcat ("Blank input, using default value: ",task_ID,"\n"));
else
  task_ID = char(prpt_taskID);
endif

%iRT data filename
prpt_iRT_fname = inputdlg ("Insert iRT unpackaged data filename. Ex: iRT data.xlsx", "Input iRT unpackaged data filename");
if (isempty ( char (prpt_iRT_fname) ) == 1 )
  printf ( cstrcat ("Blank input, using default value: ",filename,"\n"));
else
  filename = char(prpt_iRT_fname);
endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIRST PART: READING FILE

% read data from .xlsx
printf("Reading file...\n")
xls = xlsopen(filename);
session_arr = xls2oct(xls,"sessions");
raw_arr = xls2oct(xls,"data");
[xls] = xlsclose(xls);

printf("Reading data...\n");

% get reference data

for n = 1 :size(session_arr,1)
  if and( strcmp(session_arr{n,1},num2str(session_ID)), strcmp(session_arr{n,2}, task_ID) )
    ref_quat_data = cell2mat(session_arr(n,7:10));
    printf("Reference orientation found: %d %d %d %d\n",ref_quat_data(1),ref_quat_data(2),ref_quat_data(3),ref_quat_data(4));
    break;
  endif
endfor

% find slice of data corresponding to session_ID and task_ID
% find first line of the slice
first_line= -1;
for n = 1 :size(raw_arr,1)
  if and( strcmp(raw_arr{n,1},num2str(session_ID)), strcmp(raw_arr{n,2}, task_ID) )
    first_line = n;
    printf("First line of slice is %i\n", first_line);
    break;
  endif
endfor
% find last line of the slice
last_line= -1;
for n = first_line : size(raw_arr,1)
  if not(and( strcmp(raw_arr{n,1},num2str(session_ID)), strcmp(raw_arr{n,2}, task_ID) ))
    last_line = n-1;
    printf("Last line of slice is %i\n", last_line);
    break;
  endif
  if n == size(raw_arr,1)
    last_line = n;
    printf("Last line of slice is %i\n", last_line);
  endif
endfor

% get slice of data
matrix_quat_data=cell2mat(raw_arr(first_line:last_line,5:8));
% Size of the series:
N = last_line - first_line + 1;
printf("Size: %i\n", N);
% Cleaning everything except the variables...
clear -x N matrix_quat_data ref_quat_data;



% CONCLUSION: quaternions are in the array 'matrix_quat_data' (N x 4)
% Each line has four numbers, the four quaternionic coordinates
% The last one is the real part
% The quaternion of the reference object is in 'ref_quat_data' (1x4)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SECOND PART: normalizing data
% The vector part (first 3 coordinates) points to the rotation axis
% It's norm is the sine of theta/2 (theta is the rotation, between 0 and 180)
% We change the norm to be proportional to theta/180
% The 3 new coordinates will be kept inside the unidimensional arrays x, y and z

printf("Normalizing data...\n")


% Every quaternion should be an unitary quaternion, but the norm fluctuates
% Normalizing the reference quaternion
ref_quat_data = ref_quat_data / norm(ref_quat_data);
% Setting the matrix that will place reference as the identity (the 'zero')
a_ref = ref_quat_data(1);
b_ref = ref_quat_data(2);
c_ref = ref_quat_data(3);
d_ref = ref_quat_data(4);
ref_matrix_transposed = [d_ref, c_ref, -b_ref, a_ref; -c_ref, d_ref, a_ref, b_ref; b_ref, -a_ref, d_ref, c_ref; -a_ref, -b_ref, -c_ref, d_ref];

% Creating x, y, z
x = zeros(N);
y = zeros(N);
z = zeros(N);

for i = 1:N
  % The first step is normalize the quaternions
  q = matrix_quat_data(i,1:4);
  q_norm = norm(q);
  q = q / q_norm;
  % Placing reference position as the identity (the 'zero')
  q = q * ref_matrix_transposed;
  % Sign correction: real part must be positive
  if q(4) < 0
    q = -q;
  endif
  % Extracting angle of rotation (angle distance to reference) from the real part
  % Normalized to [0,1]
  angle = 360 * acos(q(4)) / pi; % angle for resolugram
  % Correction factor: the angle distance divided by the norm of the vector part
  correction_factor = angle / norm( q(1:3) );
  % Defining x , y, z
  x(i) = correction_factor * q(1);
  y(i) = correction_factor * q(2);
  z(i) = correction_factor * q(3);
endfor

% Conclusion: x, y, z are the coordinates of the worm-path inside the quaternionic ball

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PART 3: Auxiliary functions
% These functions are used to draw a path with thickness, in order to give a sense of perspective


function v_most = most_distant_vector(v,VX,VY,VZ)
  % v is a (N,3) array with the tangent vectors.
  % These vectors must be unitary or null.
  % VX,VY,VZ are arrays (2d) of unit vectors (typically, generated by
  % a grid in spherical coordinates)
  % The function chooses the unit vector which is one of the farthest from
  % the normalized vectors of the list v (except those which are zero)
  n = size(VX)(1);
  m = size(VX)(2);
  bigger = 0;
  i_final = 1;
  j_final = 1;
  for i=1:n
    for j=1:m
      dist2 = min( (v(:,1) - VX(i,j)).^2 + (v(:,2) - VY(i,j)).^2 + (v(:,3) - VZ(i,j)).^2 );
      if dist2 > bigger
        i_final = i;
        j_final = j;
        bigger = dist2;
      endif
    endfor
  endfor
  v_most = [ VX(i_final,j_final) VY(i_final,j_final) VZ(i_final,j_final)];
endfunction

function A = annotate(x,y,z,d)
  % d is the threshold to decide when an antipodal transition happens
  % x, y, z are the arrays with the coordinates of the trajectory
  % A will be a vector with integer entries 0, 1, 2, 3, 4, 5, 6, 7 or 8
  % according to the following:
  N = size(x)(1);
  A = zeros(N,1);
  for i=1:N
    %i
    go = false;
    if i == 1
      go = true;
    else
      dist_previous = norm( [ x(i)-x(i-1), y(i)-y(i-1), z(i)-z(i-1) ]);
      if (dist_previous > d)
        go = true;
      endif
    endif
    stop = false;
    if i == N
      stop = true;
    else
      dist_next = norm( [ x(i)-x(i+1), y(i)-y(i+1), z(i)-z(i+1) ]);
      if (dist_next > d)
        stop = true;
      endif
    endif
    null_dist_previous = false;
    null_dist_next = false;
    if !go
      null_dist_previous = (dist_previous == 0);
    endif
    if !stop
      null_dist_next = (dist_next == 0);
    endif
    conditions = [ go , stop , null_dist_previous , null_dist_next];
    if (conditions == [ false, false, false, false ])
      A(i) = 0;
    endif
    if (conditions == [ true, false, false, false ])
      A(i) = 1;
    endif
    if (conditions == [ false, true, false, false ])
      A(i) = 2;
    endif
    if (conditions == [ true, true, false, false ])
      A(i) = 3;
    endif
    if (conditions == [ false, false, false, true ])
      A(i) = 4;
    endif
    if (conditions == [ false, false, true, false ])
      A(i) = 5;
    endif
    if (conditions == [ false, false, true , true ])
      A(i) = 6;
    endif
    if (conditions == [ true, false, false, true ])
      A(i) = 7;
    endif
    if (conditions == [ false, true, true, false ])
      A(i) = 8;
    endif
  endfor
endfunction

function tangent_out = calculate_tangent_vectors(x,y,z,A)
  % x,y,z are the arrays with the coordinates of the trajectory
  % A is the array of the same size which is annotaded as in the function annotate
  % v will be an (N,3) array (where N is the size of x,y,z and A)
  N = size(x)(1);
  tangent_out = zeros(N,3);
  % Populating tangent_out: tangent direction departing from each point
  % It will be zero at the end of each snippet
  for i=1:N
    if ismember(A(i),[2,3,4,6,7,8])
      tangent_out(i,:) = [ 0 , 0 , 0 ];
    else
      difference_vector = [ x(i+1)-x(i), y(i+1)-y(i), z(i+1)-z(i) ];
      tangent_out(i,:) = difference_vector/norm(difference_vector);
    endif
  endfor
endfunction



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 4: Drawing axes, meridians, equator etc

printf("Drawing axes, meridians, equator...\n");

lsrc = [0 0 10];
% lsrc (light source): either a 2-element vector [azimuth, elevation]
% or a 3-element vector [lx,ly,lz]
P = [0.7 0.6 0.7 10];
% P = [AM D SP exp]
%"AM" strength of ambient light
%"D" strength of diffuse reflection
%"SP" strength of specular reflection
%"EXP" specular exponent
% Circles, as tori



axis("equal")
axis([-180 180 -180 180 -180 180])
set(gca, 'xtick', -180:90:180);
set(gca, 'ytick', -180:90:180);
set(gca, 'ztick', -180:90:180);
set(gca, 'xticklabel',({'-180','','0','','+180'}));
set(gca, 'yticklabel',({'-180','','0','','+180'}));
set(gca, 'zticklabel',({'-180','','0','','+180'}));
% perspective view
view([360 -180 250]);
% view of xy projection
%view([0 0 400]);
% view of xz projection
%view([0 -400 0]);
% view of yz projection
%view([400 0 0]);
colormap(copper(64));
%axis("Projection","perspective")

xlabel("x");
ylabel("y");
zlabel("z");

hold on

% Creation of a torus which will be used for meridians, equator and latitude lines
r = 0.003;
R = 1.0;
phi = linspace(0,2*pi,51);
theta = linspace(0,2*pi,51);
[PHI,THETA] = meshgrid(phi,theta);
% Thin torus in the xy-plane: equator
TX = 180 * (R + r*cos(THETA)).*cos(PHI);
TY = 180 * (R + r*cos(THETA)).*sin(PHI);
TZ = 180 * r*sin(THETA);

% Rotating original around x-axis
% Meridian in the xz-plane
M = [1 0 0; 0 0 1; 0 -1 0];
TXX = M(1,1)*TX + M(1,2)*TY + M(1,3)*TZ;
TYY = M(2,1)*TX + M(2,2)*TY + M(2,3)*TZ;
TZZ = M(3,1)*TX + M(3,2)*TY + M(3,3)*TZ;
%surfl(TXX,TYY,TZZ,lsrc,P);

% Meridians
angle_fraction = 2;
for i = 1:angle_fraction-1
  angle = i*(pi/angle_fraction);
  M = [ cos(angle) sin(angle) 0; -sin(angle) cos(angle) 0; 0 0 1];
  TXXX = M(1,1)*TXX + M(1,2)*TYY + M(1,3)*TZZ;
  TYYY = M(2,1)*TXX + M(2,2)*TYY + M(2,3)*TZZ;
  TZZZ = M(3,1)*TXX + M(3,2)*TYY + M(3,3)*TZZ;
  surfl(TXXX,TYYY,TZZZ,lsrc,P);
endfor

% Latitude lines
angle_fraction = 2;
for i = 1:angle_fraction-1
  angle = -pi/2 + i * (pi/angle_fraction);
  Radius = R*cos(angle)
  TXXX = 180 * (Radius + r*cos(THETA)).*cos(PHI);
  TYYY = 180 * (Radius + r*cos(THETA)).*sin(PHI);
  TZZZ = 180* ( R*sin(angle) + r*sin(THETA));
  surfl(TXXX,TYYY,TZZZ,lsrc,P);
endfor


% AXES
[CX,CY,CZ]=cylinder([1 1],16);
% axis z
M = [1 0 0; 0 1 0; 0 0 1];
CXX = 180* ( M(1,1)*r*CX + M(1,2)*r*CY + M(1,3)*(-1+2*CZ) );
CYY = 180* ( M(2,1)*r*CX + M(2,2)*r*CY + M(2,3)*(-1+2*CZ) );
CZZ = 180* ( M(3,1)*r*CX + M(3,2)*r*CY + M(3,3)*(-1+2*CZ) );
surfl(CXX,CYY,CZZ,lsrc,P);
% axis y
M = [1 0 0; 0 0 1; 0 -1 0];
CXX = 180 * ( M(1,1)*r*CX + M(1,2)*r*CY + M(1,3)*(-1+2*CZ));
CYY = 180 * ( M(2,1)*r*CX + M(2,2)*r*CY + M(2,3)*(-1+2*CZ));
CZZ = 180 * ( M(3,1)*r*CX + M(3,2)*r*CY + M(3,3)*(-1+2*CZ));
surfl(CXX,CYY,CZZ,lsrc,P);
% axis x
M = [0 0 1; 0 1 0; -1 0 0];
CXX = 180 * ( M(1,1)*r*CX + M(1,2)*r*CY + M(1,3)*(-1+2*CZ));
CYY = 180 * ( M(2,1)*r*CX + M(2,2)*r*CY + M(2,3)*(-1+2*CZ));
CZZ = 180 * ( M(3,1)*r*CX + M(3,2)*r*CY + M(3,3)*(-1+2*CZ));
surfl(CXX,CYY,CZZ,lsrc,P);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PART 5: DRAWING THE TRAJECTORY

% This is interesting to better see the antipodal transitions
%plot3(x,y,z);




printf("Drawing trajectory:\n");

% Unit sphere
phi = linspace(0,2*pi,9);
theta = linspace(-pi/2,pi/2,9);
[PHI,THETA] = meshgrid(phi,theta);
SX = cos(THETA).*cos(PHI);
SY = cos(THETA).*sin(PHI);
SZ = sin(THETA);

% radius_0 is the radius of the trajectory at the beginning
% radius_1 at the end
radius_0 = 2.0;
radius_1 = 0.1;




% Calculating tangent vectors
% Threshold to determine an antipodal transition:
d = 300
% Annotate transitions and freezings along the trajectory
printf("Parsing trajectory...\n")
A = annotate(x,y,z,d);
% Calculating tangent vectors at each point of the trajectory
printf("Tangent vectors...\n")
v = calculate_tangent_vectors(x,y,z,A);
% Obtaining the most distant unitary vector to all tangent vectors
[VX,VY,VZ] = sphere(25);
v_most = most_distant_vector(v,VX,VY,VZ);
%v_most

printf("Calculating and plotting cylinders...\n")



% First and second normal vectors 1) fn = v_most x v
% 2) sn = v x fn
fn = zeros(N,3);
sn = zeros(N,3);

for i = 1:N
  vi = [v(i,1),v(i,2),v(i,3)];
  if vi != [0,0,0]
    fni = cross(v_most,vi);
    vi = vi/norm(vi);
    fni = fni/norm(fni);
    sni = cross(vi,fni);
    fn(i,:) = fni;
    sn(i,:) = sni;
  endif
endfor

counter = 0
for i=1:N-1
  this_t = sqrt(i/N);
  next_t = sqrt((i+1)/N);
  % radius is the linear interpolation between radius_0 and radius_1
  this_radius = (1-this_t)*radius_0 + this_t*radius_1;
  next_radius = (1-next_t)*radius_0 + next_t*radius_1;
  radii = [this_radius,next_radius]';
  surfl(x(i)+this_radius*SX,y(i)+this_radius*SY,z(i)+this_radius*SZ,lsrc,P);
  if ismember(A(i),[1,2,3,7,8])
    surfl(x(i)+2*this_radius*SX,y(i)+2*this_radius*SY,z(i)+2*this_radius*SZ,lsrc,P);
    label_number = num2str(counter);
    text(x(i)+2*this_radius,y(i)+2*this_radius,z(i)+4*this_radius, label_number, 'fontsize', 20);
    counter = counter + 1;
  endif
  if ismember(A(i),[0,1,5])
    CXX = [x(i),x(i+1)]' + radii.*[fn(i,1),fn(i,1)]'.*CX + radii.*[sn(i,1),sn(i,1)]'.*CY;
    CYY = [y(i),y(i+1)]' + radii.*[fn(i,2),fn(i,2)]'.*CX + radii.*[sn(i,2),sn(i,2)]'.*CY;
    CZZ = [z(i),z(i+1)]' + radii.*[fn(i,3),fn(i,3)]'.*CX + radii.*[sn(i,3),sn(i,3)]'.*CY;
    surfl(CXX, CYY, CZZ, lsrc,P);
    shading interp;
    hold on
  endif
endfor

printf("\nProcess completed!\n"); toc();

