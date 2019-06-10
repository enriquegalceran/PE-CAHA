import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from Salida_limpia import mostrarresultados, stdrobusta
from astropy.io import fits
import IMGPlot as ImP
import json


def deshacer_tupla_coord(tupla):
    return tupla[0], tupla[1], tupla[2], tupla[3]


# thisdict =	{
#   "brand": (1, "Ford"),
#   "model": (2, "Mustang"),
#   "year": (3, 1964)
# }
# print(thisdict)


# sacar keywords
# for x in thisdict:
#   print(x)

# sacar valores
# for x in thisdict:
#   print(thisdict[x])

# sacar valores alt
# for x in thisdict.values():
#   print(x)

##########################################################################
# El objetivo es crear una base de datos (diccionario) con los nombres de los filtros

def mostrar_diccionario(nombre_diccionario):
    print("Mostramos el diccionario: ")
    for x, y in nombre_diccionario.items():
      print(x, y)


def guardar_json(variable, filename):
    json_f = json.dumps(variable)
    f = open(filename,"w")
    f.write(json_f)
    f.close()


def cargar_json(filename='Dic_filtro.json'):
    with open(filename) as json_file:
        data = json.load(json_file)
    return data


Dic_filtro = cargar_json()
Dic_filtro = {
  "free": 0
}

nombre = "CousinsR"


def leer_diccionario(nombre, diccionario_=Dic_filtro, filename='Dic_filtro.json'):
    if nombre in diccionario_:
        return diccionario_[nombre]
    else:
        len_dic = len(diccionario_)
        diccionario_[nombre] = len_dic
        print('Nueva entrada en el diccionario: ', diccionario_[nombre])
        guardar_json(diccionario_, filename)
        return diccionario_[nombre]


valor = leer_diccionario("CousinsR")
print(valor)
print(type(Dic_filtro)==dict)
mostrar_diccionario(Dic_filtro)

# prueba = str(Dic_filtro['free'])
# print(prueba)
# print(len(prueba))

