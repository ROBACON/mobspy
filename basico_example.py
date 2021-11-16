import basico
from models import builder
from tools import rewrite_cps
from utilities import parse
import logging
import run
import tempfile
import copy
import os
import shutil
import codecs

logging.basicConfig()
module_logger = logging.getLogger('root')

params = parse.params("AB_2CD.txt", module_logger)

temp = tempfile.mkdtemp()
temp += '/'

instance_params = copy.deepcopy(params)
instance_params['cps_dir'] = temp
instance_params['output_dir'] = temp
instance_params['output_file'] = os.path.join(temp, os.path.basename(params['output_file']))
instance_params['save_data'] = True

mappings = run.init(instance_params)
print(instance_params['cps_dir'])
print(os.listdir(instance_params['cps_dir']))
exit()
data = run.run_solver(instance_params, mappings)
sbml_file = instance_params['cps_dir'] + "/model.sbml"

sbml_str = ""
with codecs.open(sbml_file, "r", "utf-8") as f:
    for line in f.readlines():
        sbml_str += line

model = basico.model_io.load_model_from_string(sbml_str)
ans = basico.run_time_course()
print(ans)

# print(data)
shutil.rmtree(temp)