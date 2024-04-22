import json
from types import SimpleNamespace
import os 

def load_paths():
    # Opening JSON file
    f = open('paths.json')
    
    # returns JSON object as 
    # a dictionary
    data = json.load(f,object_hook=lambda d: SimpleNamespace(**d))
    return data

def get_input_file(baseDir):
    fname = ""
    raw_name = ""
    counter = 0
    for file in os.listdir(baseDir):
        if file.endswith(".json"):
            fname = os.path.join(baseDir, file)
            raw_name = file
            counter+=1
    if fname == "":
        return False, fname, "Nothing to do.",raw_name
    else:
        return True, fname, "Running Model...",raw_name

def open_file(fname):
    # Opening JSON file
    f = open(fname)
    
    # returns JSON object as 
    # a dictionary
    data = json.load(f,object_hook=lambda d: SimpleNamespace(**d))
    return data