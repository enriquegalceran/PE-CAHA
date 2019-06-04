import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from Salida_limpia import mostrarresultados, stdrobusta
from astropy.io import fits
import IMGPlot as ImP
import json


def deshacer_tupla_coord(tupla):
    return tupla[0], tupla[1], tupla[2], tupla[3]


thisdict =	{
  "brand": (1, "Ford"),
  "model": (2, "Mustang"),
  "year": (3, 1964)
}
print(thisdict)


# sacar keywords
# for x in thisdict:
#   print(x)

# sacar valores
# for x in thisdict:
#   print(thisdict[x])

# sacar valores alt
# for x in thisdict.values():
#   print(x)

# Loop de ambos
# for x, y in thisdict.items():
#   print(x, y)



##########################################################################
# El objetivo es crear una base de datos (diccionario) con los nombres de los filtros

def guardar_json(variable, filename):
    json_f = json.dumps(variable)
    f = open(filename,"w")
    f.write(json_f)
    f.close()


def cargar_json(filename):
    with open('dict.json') as json_file:
        data = json.load(json_file)
    return data


guardar_json(thisdict, 'dict.json')
diccionario = cargar_json('dict.json')
print(diccionario)

if 'model' not in diccionario:
    print('no esta')
else:
    print('si esta')

print(diccionario["model"][0])