
if ispc()
vars
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\spot-master'))
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\exportfig'));
rmpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\CS_simulations\')) 
addpath(genpath('L:\basic\divi\Projects\cosart\CS_simulations\Code\general_recon'))
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\nifti\'))
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\Wavelab850\'))

else
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/spot-master'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/nifti'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/exportfig'));
addpath(genpath('/opt/amc/bart/')); vars
rmpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations')) 
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Code/general_recon'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/Wavelab850'));
addpath(genpath('/opt/amc/matlab/toolbox/MRecon/'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/imagine'));


end