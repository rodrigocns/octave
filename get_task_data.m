#returns data related to specific task number. Change or add values accordingly.
function [model_name,task_name,ref_quat] = get_task_data (task_num)
  switch task_num
  case 1
    model_name = 'pseudobatracotoxin_molecule'; #caprolactama, pseudobatracotoxin_molecule
    task_name = 'bolaBastao_c';
    ref_quat = [0.687341, -0.490258, -0.282534, -0.455395]; % quaternios da referencia
  case 2
    model_name = 'pseudobatracotoxin_molecule';
    task_name = 'poligonFill';
    ref_quat = [0.687341, -0.490258, -0.282534, -0.455395];
  case 3
    model_name = 'MRT_VK_mol';
    task_name = 'mrt';
    ref_quat = [0.921868, 0.299222, 0.223654, -0.102977];
  otherwise
    error ("ERRO! Task not recognized! Update data from get_task_data.m");
  endswitch
endfunction

