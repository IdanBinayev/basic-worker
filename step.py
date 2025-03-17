import json
import sys, traceback
import os
# data = sys.argv[1]
tr = sys.argv[1]
sys.path.insert(1, tr)
global data
with open(tr + '/get_data.json') as json_file:
  data = json.load(json_file)

def uploaded_attachment_from_trigger():
  return '../files/' + variable('attachment_from_trigger')['filename']

def log(*e):
  print(e)

def print_data():
  print(json.dumps(data, indent=4, sort_keys=True))

def variable(k):
  return data['inputs'][k]

def safe_variable(k):
  return data['inputs'].get(k)

def set_file(k, v):
  data['inputs'][k] = tr + '/rendered_files/' + v

def set_file_for_email(k, v):
  parts = tr.split('/')

  if len(parts) > 1 and parts[1] == 'tmp': # Transform the path by replacing 'tmp' with the shared path
    transformed_path = os.path.join(SHARED_TMP_PATH, *parts[2:])
  elif len(parts) > 1 and parts[1] == 'srv': # For production environments, if the path starts with 'srv', skip the first part.
    transformed_path = os.path.join(SHARED_TMP_PATH, *parts[1:])
  else:
    transformed_path = tr

  data['inputs'][k] = {'tmp_path' : os.path.join(transformed_path, 'rendered_files', v)}
  return v

def store_variable(k, v):
  data['inputs'][k] = v
  return v

def base():
  return data['base']

def token():
  return data['token']

try:
  import code
  import lab_guru_api
  from lib.LabguruPython.labguru import Labguru
  lab = Labguru(token(), base())
  lab_guru_api.setup(base(), token())

  os.chdir(tr + '/rendered_files')
  code.execution_script(lab, lab_guru_api, json, variable, store_variable, safe_variable, token, print_data, set_file, set_file_for_email, log, base, tr, uploaded_attachment_from_trigger, data)

  print('nuzi')
  print(':::EXEC SUCCESS:::')
  with open('../return_data.json', 'w') as outfile:
    json.dump(data, outfile)

except:
  e = sys.exc_info()[0]
  print( ":::EXEC ERROR::: %s" % e )
  traceback.print_exc(file=sys.stdout)
