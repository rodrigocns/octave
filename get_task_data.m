#returns data related to specific task number. Change or add values accordingly.
function [model_name,task_name,ref_quat] = get_task_data (task_num)
  switch task_num
  case 1
    model_name = 'pseudobatracotoxin_molecule'; #caprolactama, pseudobatracotoxin_molecule
    task_name = 'cor';
    ref_quat = [0.6873,  -0.4903,  -0.2825,  -0.4554]; #quaternios para a referencia
  case 2
    model_name = 'pseudobatracotoxin_molecule';
    task_name = 'cinza';
    ref_quat = [0.6873,  -0.4903,  -0.2825,  -0.4554]; #quaternios para a referencia
  case 3
    model_name = 'caprolactama';
    task_name = 'intro';
    ref_quat = [0.5183,   0.5105,   0.4575,  -0.5105]; #quaternios para a referencia
  otherwise
    error ("ERRO! Task não reconhecida! Atualize os dados de get_task_data.m");
  endswitch
endfunction
