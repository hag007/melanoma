import shutil
import os
import json
from omics.sync_omics import sync_from_project
from utils import download_resources
import sys
sys.path.insert(0, '../')
import constants

dir_path = os.path.dirname(os.path.realpath(__file__))

def main():
    dest = constants.BASE_PROFILE
    if not os.path.exists(os.path.join(dest, "output")):
        os.makedirs(os.path.join(dest, "output"))
    sync_from_project()
    download_resources.main()


if __name__ == "__main__":
    main()