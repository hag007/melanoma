import shutil
import os
import json
import sys
sys.path.insert(0, '../')
import constants

dir_path = os.path.dirname(os.path.realpath(__file__))

def sync_from_project():
    dest = constants.BASE_PROFILE
    source = dir_path


def sync_to_project():
    dest = dir_path
    source = constants.BASE_PROFILE


def replace_folder(source, dest):
    if os.path.exists(dest):
        os.rmdir(os.path.join(dest, "dictionaries"))
        os.rmdir(os.path.join(dest, "list"))
    else:
        os.makedirs(dest)
        
    shutil.copytree(os.path.join(source,"dictionaries"), os.path.join(dest, "dictionaries"))
    shutil.copytree(os.path.join(source,"list"), os.path.join(dest, "list"))
