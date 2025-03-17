import requests
import json
import random
import string
import os
from urllib.parse import urlencode

base=""
token=""
def setup(base_url, token1):
  global token
  global base
  base = base_url
  token = token1

def api_request(api, data, method='get', data_as_query_string=False):
  url = f"{base}/api/v1/{api}?token={token}"  
  if data_as_query_string:
    url += f"&{urlencode(data)}"
  try:
    if (method == 'post'):                           
      response = requests.post(url,json=data)        
    elif (method == 'put'):
      response = requests.put(url, data=data)
    elif (method == 'delete'):
      response = requests.delete(url)
    else:
      headers = {'Content-Type': 'application/json'}
      response = requests.get(url, headers=headers)
    response.raise_for_status()
    return response.json()
  except Exception as error:  
    print(error)

def get_api_data_as_query_string(api, data = {}):
  print(f"get_api_data_as_query_string({api}, {data})")
  return api_request(api, data, data_as_query_string = True)

# GET api query to laburu, parses the response and returns a hash (json-like object)
def get_api(api, data = {}):
  return api_request(api, data)

# POST api query to laburu, parses the response and returns a hash (json-like object)
def post_api(api, data):
  if data is None:
    raise ValueError("post data can't be blank")
  return api_request(api, data, method='post')

# PUT api query to laburu, parses the response and returns a hash (json-like object)
def put_api(api, data):
  if data is None:
    raise ValueError("put data can't be blank")
  return api_request(api, data, method='put')

# DELETE api query to laburu, parses the response and returns a hash (json-like object)
def delete_api(api):
  return api_request(api, {}, method='delete')

def link_objects(source_uuid, target_uuid):
  print(f"link_objects({source_uuid}, {target_uuid}")
  return post_api('links.json', {"item": {"source_uuid": source_uuid, "target_uuid": target_uuid}})

def create_experiment_from_protocol(title, project_id, protocol_id, milestone_id):
  data = {'item': {'title': title, 'project_id': project_id, 'milestone_id': milestone_id, 'protocol_id': protocol_id}}
  print(f"create_experiment_from_protocol({title}, {project_id},{protocol_id}, {milestone_id})")
  return post_api("experiments.json",data)

def get_user(id):
  print(f"get_user({id})")
  return get_api(f"users/{id}.json")

def get_all_protocols():
  print(f"get_all_protocols")
  return get_api("protocols.json")

def get_all_folders(project_id = {}):
  print(f"get_all_folders({project_id})")
  ret = get_api("milestones.json")
  if not project_id:
    return ret
  else:
    response = []
    for folder in ret:
      if folder['project_id'] == project_id:
        response.append(folder)
    return response

def create_folder(project_id, title, description):
  print(f"create_folder({project_id}, {title}, {description})")
  data = {'item' : { 'title' : title, 'description' : description, 'project_id' : project_id }}
  return post_api("milestones.json", data)


def get_all_attachments():
  try:
    # &filter="+"{"+ f"attachable_type:Projects::Experiment,attachable_id:{id}"+"}
    print(f"get_all_attachments()")
    return get_api(f"/attachments")
  except:
    return "" 

def get_all_projects():
  print('get_all_projects()')
  return get_api('projects.json')
                        
# query project by id
def get_project(id):
  if id is None:
    raise ValueError("Project#id can't be blank")
  print(f"get_project({id})")
  return get_api(f'projects/{id}.json')

 # query protocol by id
def get_protocol(id):  
  if id is None:
    raise ValueError("Protocol#id can't be blank")
  print(f"get_protocol({id})")
  return get_api(f"protocols/{id}.json")

# query section by id
def get_section(id):
  if id is None:
    raise ValueError("Section#id can't be blank")
  print(f"get_section({id})")
  return get_api(f"sections/{id}.json")

# query element by id
def get_element(id):
  if id is None:
    raise ValueError("Element#id can't be blank")
  print(f"get_element({id})")
  return get_api(f"elements/{id}.json")

# query report by id
def get_report(id):
  if id is None:
    raise ValueError("Report#id can't be blank")
  print(f"get_report({id})")
  res = get_api(f"reports")
  response = []
  for report in res:
    if report['id'] == id:
      return report

# query report by name
def get_report_by_name(name):
  if name is None:
    raise ValueError("Report#name can't be blank")
  res = get_api("reports.json")
  response = []
  for report in res:
    if report['title'] == name:
      return report

def big_rand():
  source = string.ascii_letters + string.digits
  return ''.join((random.choice(source) for i in range(36)))

def tmp_file_generator(ext):
  return f"./tmp/{big_rand()}" + (('.' + ext) if ext else '')

# downloads an attachment from labguru by attachment id and file [ext]ension
def download_attachment(*args):
  id, ext = args[0]['id'], args[0]['extension']
  ret = requests.get(f"{base}/api/v1/attachments/{id}/download?token={token}", allow_redirects=True)  
  f = open('../../.'+tmp_file_generator(ext), 'wb')
  f.write(ret.content) 
  f.close()
  print(f"download_attachment({id}, {ext})")
  return f

# uploads an attachment to labguru
def upload_attachment(title, raw_attachment, description):
  print(f"upload_attachment({title}, {raw_attachment}, {description})")
  return post_api("attachments.json", { "item": { "title": title, "attachment": raw_attachment, "description": description }})

# query attachment by id (attachment meta data, not the raw attachment)
def get_attachment(id):
  if id is None:
    raise ValueError("Attachment#id can't be blank")
  print(f"get_attachment({id})")
  return get_api(f"attachments/{id}.json")

def get_document(id):
  if id is None:
    raise ValueError("Document#id can't be blank") 
  print(f"get_document({id})")
  return get_api(f"documents/{id}")

# query experiment by id
def get_experiment(id):
  if id is None:
    raise ValueError("Experiment#id can't be blank")
  print(f"get_experiment({id})")
  return get_api(f"experiments/{id}.json")

# create a report with a title
def create_report(title):
  if title is None:
    raise ValueError("Report#title can't be blank")
  print(f"create_report({title})")  
  return post_api("reports.json", {"item": {"title": title }})

# create an element with content as [data] and container_id
def create_element(data, container_id):
  print(f"create_element({data}, {container_id})")
  return post_api("elements.json",{ "item": { "element_type": "text",
                                              "data": data,
                                              "container_id": container_id,
                                              "container_type": "ExperimentProcedure"}})