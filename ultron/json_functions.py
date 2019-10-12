import json


def load_json(filename='Dic_filtro.json'):
    with open(filename) as json_file:
        data = json.load(json_file)
    return data


def save_json(variable, filename):
    json_f = json.dumps(variable, indent=2, sort_keys=True)
    f = open(filename, "w")
    f.write(json_f)
    f.close()
