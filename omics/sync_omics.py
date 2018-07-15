import shutil
import os
import json

dir_path = os.path.dirname(os.path.realpath(__file__))
config_json = json.load("config/conf.json")

def sync_from_project():
    dest = config_json['BASE_PROFILE']
    source = dir_path


def sync_to_project():
    dest = dir_path
    source = config_json['BASE_PROFILE']


def replace_folder(source, dest):
    if os.path.exists(dest):
        os.rmdir(os.path.join(dest, "dictionaries"))
        os.rmdir(os.path.join(dest, "list"))
    else:
        os.makedirs(dest)
        
    shutil.copytree(os.path.join(source,"dictionaries"), os.path.join(dest, "dictionaries"))
    shutil.copytree(os.path.join(source,"list"), os.path.join(dest, "list"))
