# -*- coding: utf-8 -*-

"""
    Proyecto de Unidad de Limpieza y Tratamiento de Resultados Observacionales Nativo (Proyecto U.L.T.R.O.N.)

    Recoge los resultados de observaciones del instrumento CAFOS y los reduce correctamente.
    Desarrollado por Enrique Galceran García
"""

from astropy.io import fits
from Salida_limpia import mostrarresultados, stdrobusta
import numpy as np
import numpy.ma as ma
import IMGPlot as ImP
# import matplotlib.pyplot as plt
import os
import csv
import argparse
import time
import warnings
import json


def deshacer_tupla_coord(tupla):
    return tupla[0], tupla[1], tupla[2], tupla[3]


def guardar_listas_csv(csvfile, res):
    with open(csvfile, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        for val in res:
            writer.writerow([val])


def mostrar_diccionario(nombre_diccionario):
    print("Mostramos el diccionario: ")
    for x, y in nombre_diccionario.items():
        print(x, y)


def guardar_json(variable, filename):
    json_f = json.dumps(variable, indent=2, sort_keys=True)
    f = open(filename, "w")
    f.write(json_f)
    f.close()


def cargar_json(filename='Dic_filtro.json'):
    with open(filename) as json_file:
        data = json.load(json_file)
    return data


def conseguir_listas_archivos(path_):
    listaarchivos = []
    for file in os.listdir(path_):
        if file.endswith(".fits"):
            listaarchivos.append(os.path.join(path_, file))
    return listaarchivos


def leer_diccionario(nombre, diccionario_, filename='Dic_filtro.json'):
    """
    Busca en el diccionario el número del filtro buscado. Si no existe, crea uno nuevo

    :param nombre:
    :param diccionario_:
    :param filename:
    :return:
    """
    if nombre in diccionario_:
        return diccionario_[nombre]
    else:
        len_dic = len(diccionario_)
        diccionario_[nombre] = len_dic
        print('Nueva entrada en el diccionario: Filtro {0} - Indice {1:03d}'.format(nombre, diccionario_[nombre]))
        guardar_json(diccionario_, filename)
        return diccionario_[nombre]


def sacar_coordenadas_ccd(imagen_, mypath_=False):
    """
    Dado un string de texto en formato de CAFOS nos devuelve los valores x1, x2, y1 e y2.

    :param imagen_:
    :param mypath_:
    :return:
    """
    if mypath_:
        str_ccdsec = fits.open(mypath_ + imagen_)[0].header['CCDSEC']
        longitud = len(str_ccdsec)
    else:
        str_ccdsec = imagen_
        longitud = len(str_ccdsec)+1

    coma = None
    for x in range(1, longitud):
        if str_ccdsec[x] == ',':
            coma = x
            break
    if coma is None:
        raise ValueError('coma not defined!')

    puntos = None
    for x in range(coma + 1, longitud):
        if str_ccdsec[x] == ':':
            puntos = x
            break
    if puntos is None:
        raise ValueError('puntos not defined!')

    coma2 = None
    for x in range(puntos + 1, longitud):
        if str_ccdsec[x] == ',':
            coma2 = x
            break
    if coma2 is None:
        raise ValueError('coma2 not defined!')

    x1 = int(str_ccdsec[1:coma])
    y1 = int(str_ccdsec[coma + 1: puntos])
    x2 = int(str_ccdsec[puntos + 1: coma2])
    y2 = int(str_ccdsec[coma2 + 1: -1])
    tupla_salida = (x1, x2, y1, y2)
    return tupla_salida


def leer_lista(archivo):
    with open(archivo, 'rt') as f:
        reader = csv.reader(f, delimiter=',')
        your_list = list(reader)
    your_list = [item for sublist in your_list for item in sublist]
    return your_list


def sacar_coordenadas_2(lista, idx):
    x1 = lista[idx, 0]
    x2 = lista[idx, 1]
    y1 = lista[idx, 2]
    y2 = lista[idx, 3]
    tupla_salida = (x1, x2, y1, y2)
    return tupla_salida


def obtener_bias(dir_bias_, noche, lista_noches, lista_bias, x1, x2, y1, y2, b1, b2,
                 busqueda_max=10, relleno=680, verbose=False):

    """
    Busca en la lista de bias disponibles si exite un archivo de bias para la misma noche en la que se hizo la foto.
    Si exite, usará esa directamente. Si no existe, busca una imagen de bias diferente del mismo autor de noches
    cercanas. Si no hay, genera un bias plano.

    :param dir_bias_:
    :param noche:
    :param lista_noches:
    :param lista_bias:
    :param x1:
    :param x2:
    :param y1:
    :param y2:
    :param b1:
    :param b2:
    :param busqueda_max:
    :param relleno:
    :param verbose:
    :return:
    """

    bias_asociado_nombre = dir_bias_ + noche +\
        "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}.fits".format(x1, x2, y1, y2, b1, b2)
    existe = None
    cogido_otro = False
    if bias_asociado_nombre in lista_bias:
        existe = True
    else:
        posicion = lista_noches.index(noche)
        for i in range(1, busqueda_max):
            for mult in [-1, 1]:
                indice = i * mult
                noche = lista_noches[posicion + indice]
                bias_asociado_nuevo = dir_bias_ + noche + \
                    "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}.fits"\
                    .format(x1, x2, y1, y2, b1, b2)

                if bias_asociado_nuevo in lista_bias:
                    bias_asociado_nombre = bias_asociado_nuevo
                    existe = True
                    cogido_otro = True

    if existe:
        if verbose:
            if cogido_otro:
                print('Se ha cogido un bias de un dia diferente. El bias que se ha cogido es:')
                print(bias_asociado_nombre)
            else:
                print('Existe el bias asociado')
        bias_asociado = fits.getdata(bias_asociado_nombre, ext=0)
    else:
        if verbose:
            print('No hay bias cercanos validos. Se ha generado uno.')
        naxis1_expected = int((y2 - y1 + 1) / b2)
        naxis2_expected = int((x2 - x1 + 1) / b1)
        bias_asociado = np.full((naxis1_expected, naxis2_expected), relleno, dtype=float)

    return existe, bias_asociado_nombre, bias_asociado


def create_circular_mask(h, w, center=None, radius=None):

    """
    Genera la mascara circular en función de la altura y anchura de la imagen además del centro y radio.

    :param h:
    :param w:
    :param center:
    :param radius:
    :return:
    """

    if center is None:  # use the middle of the image
        center = [int(w/2), int(h/2)]
    if radius is None:  # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    yc, xc = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((xc - center[0])**2 + (yc-center[1])**2)

    mask = dist_from_center <= radius
    return mask


def crear_lista_unicos(dir_datos, noche, lista_cosas, cabecera, binning=False, nombre_filtro=False):

    """
    Dado una lista de parámetros, genera una serie de listas:
        - Lista de archivos
        - Lista de valores únicos (sin repetidos)
        - Para cada valor único, las instancias que se repiten
        - Para cada elemento de la lista, genera un índice de a cual de ĺos valores únicos pertenece
        - Si se pide, la cantidad de secciones que tiene
        - Si se pide, el nombre de los filtros en cuestión

    :param dir_datos:
    :param noche:
    :param lista_cosas:
    :param cabecera:
    :param binning:
    :param nombre_filtro:
    :return:
    """
    lista = []

    for imagen in lista_cosas:
        lista.append(fits.open(dir_datos + noche + '/' + imagen)[0].header[cabecera])
    lista_unicas, lista_count = np.unique(lista, return_counts=True)  # Contamos secciones unicas

    bin_secciones = -1 * np.ones((len(lista_unicas), 2), dtype=int)
    indice_cosas = np.zeros(len(lista_cosas), dtype=int)
    nombres_filtros = []

    for i in range(len(lista_cosas)):
        for j in range(len(lista_unicas)):
            if lista[i] == lista_unicas[j]:
                indice_cosas[i] = j
                if binning:
                    bin_secciones[j, 0] = int(fits.open(dir_datos + noche + '/' + lista_cosas[i])[0].header['crpix2'])
                    bin_secciones[j, 1] = int(fits.open(dir_datos + noche + '/' + lista_cosas[i])[0].header['crpix1'])
                if nombre_filtro:
                    nombres_filtros.append(fits.open(dir_datos + noche + '/' + lista_cosas[i])[0].header['INSFLNAM'])
                break

    return lista, lista_unicas, lista_count, indice_cosas, bin_secciones, nombres_filtros


def listas_archivos(path_):
    lista_cal = []
    lista_sci = []
    lista_misc = []
    listaarchivos = []
    for file in os.listdir(path_):
        if file.endswith(".fits"):
            if '-cal-' in file:
                lista_cal.append(file)
            elif '-sci-' in file:
                lista_sci.append(file)
            else:
                lista_misc.append(file)
            listaarchivos.append(os.path.join(path_, file))
    return lista_cal, lista_sci, lista_misc, listaarchivos


def imagen_mas_probable(archivo):
    image_data = fits.getdata(archivo, ext=0)
    v_median = np.median(image_data)
    v_sd = stdrobusta(image_data)
    if v_sd < 50 and v_median < 900:
        probable = 0  # Bias
    elif v_sd < 500:
        probable = 1  # Arc
    else:
        probable = 2  # Flat

    return probable


def comprobar(archivo, descriptor, descriptor2=None, descriptor3=None, verbose=False):
    if any([descriptor2, descriptor3]):
        descriptortemp = descriptor + descriptor2 + descriptor3
        coincide = True
        for texto in descriptortemp:
            if texto in fits.open(archivo)[0].header['OBJECT']:
                # Alguno coincide, no es ciencia
                coincide = False
                break
    else:
        coincide = False
        for texto in descriptor:
            if texto in fits.open(archivo)[0].header['OBJECT']:
                # Sale bien
                coincide = True
                break
    if verbose:
        print(archivo, fits.open(archivo)[0].header['OBJECT'], fits.open(archivo)[0].header['imagetyp'], coincide)

    # Quiza meterle algo de que si object=test que lo ignore

    return coincide


def listas_archivos2(path_, desc_bias, desc_flats, desc_arc, verbose=False, calysci=True):
    lista_cal = []
    lista_sci = []
    lista_misc = []
    lista_bias = []
    lista_flat = []
    lista_ciencia = []
    lista_arc = []
    lista_else = []
    listaarchivos = []
    lista_falla = []

    # Miramos todos los archivos dentro de la carpeta
    for file in os.listdir(path_):
        if file.endswith(".fits"):
            if calysci:
                if '-cal-' in file:
                    lista_cal.append(file)
                elif '-sci-' in file:
                    lista_sci.append(file)
                else:
                    lista_misc.append(file)
            listaarchivos.append(os.path.join(path_, file))

            # Separamos segun IMAGETYP
            tipo = fits.open(path_ + file)[0].header['IMAGETYP'].strip()
            if not fits.open(path_ + file)[0].header['OBJECT'] == 'Test':  # Comprobamos que no es un test
                if tipo == 'BIAS' or tipo == 'bias':
                    coincide = comprobar(path_ + file, desc_bias, verbose)
                    if coincide:
                        lista_bias.append(file)
                    else:
                        lista_falla.append(file)

                elif tipo == 'flat':
                    coincide = comprobar(path_ + file, desc_flats, verbose)
                    if coincide:
                        lista_flat.append(file)
                    else:
                        lista_falla.append(file)

                elif tipo == 'arc':
                    coincide = comprobar(path_ + file, desc_arc, verbose)
                    if coincide:
                        lista_arc.append(file)
                    else:
                        lista_falla.append(file)

                elif tipo == 'science':
                    coincide = comprobar(path_ + file, desc_bias, desc_flats, desc_arc, verbose)
                    if coincide:
                        lista_ciencia.append(file)
                    else:
                        lista_falla.append(file)

                else:
                    lista_else.append(file)

    for file in lista_falla:

        if calysci:
            if file in lista_sci:
                lista_ciencia.append(file)

        probable = imagen_mas_probable(path_ + file)

        if probable == 0:
            lista_bias.append(file)
        elif probable == 1:
            lista_arc.append(file)
        elif probable == 2:
            lista_flat.append(file)

    # Aseguramos que estan ordenadas
    # lista_bias = lista_bias.sort()
    # lista_flat = lista_flat.sort()
    # lista_arc = lista_arc.sort()
    # lista_ciencia = lista_ciencia.sort()
    # listaarchivos = listaarchivos.sort()
    # lista_falla = lista_falla.sort()

    return lista_bias, lista_flat, lista_arc, lista_ciencia, listaarchivos, lista_falla


def crear_listas_cal_y_sci(lista_noches_, dir_listas_, dir_datos_, desc_bias, desc_flats, desc_arc, verbose, calysci):
    i = 0
    for noche in lista_noches_:
        i += 1
        if noche not in os.listdir(dir_listas_):
            os.mkdir(dir_listas_ + noche + '/')
            # Quiza hay que descomentar para que funcione
            # CAL = []
            # SCI = []
            # MIS = []
            # ARC = []
            path_ = dir_datos_ + noche + '/'
            l_bias, l_flat, l_arc, l_ciencia, l_archivos, l_falla = listas_archivos2(path_,
                                                                                     desc_bias,
                                                                                     desc_flats,
                                                                                     desc_arc,
                                                                                     verbose,
                                                                                     calysci)

            mostrarresultados(['Bias', 'Flat', 'Arc', 'Ciencia', 'Falla'],
                              [len(l_bias), len(l_flat), len(l_arc), len(l_ciencia), len(l_falla)],
                              titulo=noche, contador=i, valor_max=len(lista_noches_))

            # guardar_listas_csv(dir_listas_ + noche + '/' + 'CAL.csv', cal)
            guardar_listas_csv(dir_listas_ + noche + '/' + 'SCI.csv', l_ciencia)
            guardar_listas_csv(dir_listas_ + noche + '/' + 'ARC.csv', l_archivos)
            guardar_listas_csv(dir_listas_ + noche + '/' + 'LBias.csv', l_bias)
            guardar_listas_csv(dir_listas_ + noche + '/' + 'LFlat.csv', l_flat)
            guardar_listas_csv(dir_listas_ + noche + '/' + 'LArc.csv', l_arc)
            guardar_listas_csv(dir_listas_ + noche + '/' + 'Errores.csv', l_falla)


def add_to_file(file, what):
    open(file, 'a+').write(what)


def slicing_data(slicing_push, size_da, size_mb):
    if slicing_push == "NW":
        s1 = size_da[0] - size_mb[0]
        s2 = size_da[0]
        s3 = 0
        s4 = size_mb[1]
    elif slicing_push == "SE":
        s1 = 0
        s2 = size_mb[0]
        s3 = size_da[1] - size_mb[1]
        s4 = size_da[1]
    elif slicing_push == "NE":
        s1 = size_da[0] - size_mb[0]
        s2 = size_da[0]
        s3 = 0
        s4 = size_mb[1]
    else:
        s1 = 0
        s2 = size_mb[0]
        s3 = 0
        s4 = size_mb[1]
    return s1, s2, s3, s4

 #
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################


def juntar_imagenes_bias(noche, secciones_unicas_, coordenadas_secciones_, secciones_count_, indice_seccion_,
                         bin_secciones_, dir_bias_, dir_datos_, lista_bias_, verbose=0,
                         interactive=False, recortar=False):
    for seccion in range(len(secciones_unicas_)):
        print('seccion: ' + str(seccion))
        coordenadas_dibujo = sacar_coordenadas_2(coordenadas_secciones_, seccion)
        x1, x2, y1, y2 = deshacer_tupla_coord(coordenadas_dibujo)

        # Sacar el Binning
        crpix1 = bin_secciones_[seccion, 1]
        crpix2 = bin_secciones_[seccion, 0]

        naxis1_expected = int((coordenadas_dibujo[3] - coordenadas_dibujo[2] + 1) / crpix1)
        naxis2_expected = int((coordenadas_dibujo[1] - coordenadas_dibujo[0] + 1) / crpix2)
        if (coordenadas_dibujo[3] - coordenadas_dibujo[2] + 1) % crpix1 != 0:
            naxis1_expected += 1
        if (coordenadas_dibujo[1] - coordenadas_dibujo[0] + 1) % crpix2 != 0:
            naxis2_expected += 1

        master_biases = np.zeros((secciones_count_[seccion],
                                  naxis1_expected,
                                  naxis2_expected), dtype=float)

        indice0 = 0
        slicing_push = False
        cabecera = None
        for imagen in range(len(lista_bias_)):
            if indice_seccion_[imagen] == seccion:
                image_file = dir_datos_ + noche + '/' + lista_bias_[imagen]
                image_data = fits.getdata(image_file, ext=0)
                if image_data[:, :].shape == master_biases[indice0, :, :].shape:
                    master_biases[indice0, :, :] = image_data[:, :]                     # Juntar
                    if indice0 == 0:
                        cabecera = fits.open(image_file)[0].header
                else:
                    size_mb = master_biases[indice0, :, :].shape
                    size_da = image_data[:, :].shape
                    if recortar:
                        if not slicing_push:
                            warnings.warn("Sizes Incompatible!")
                            print("Sizes incompatible:")
                            print("Data size: " + str(size_da))
                            print("Master Bias size: " + str(size_mb) + "\n")
                            slicing_push = (input("Slicing fits. "
                                                  "Select side towards to push (SW), SE, NE, NW: ") or "SW")
                        s1, s2, s3, s4 = slicing_data(slicing_push, size_da, size_mb)
                        master_biases[indice0, :, :] = image_data[s1:s2, s3:s4]
                    else:
                        warnings.warn("Sizes Incompatible!")
                        print("Sizes incompatible:")
                        print("Data size: " + str(size_da))
                        print("Master Bias size: " + str(size_mb) + "\n")
                        print("Skipping current Master Bias.")
                        print("Consider using slicing with '--recortar'. ")
                        input("Press Enter to continue...")
                        break
                indice0 += 1

        master_bias_colapsado = np.median(master_biases, axis=0)
        # print(master_bias_colapsado)
        # print(master_bias_colapsado.shape)
        # print(master_bias_colapsado)
        # plt.imshow(master_bias_colapsado)
        # ImP.imgdibujar(master_bias_colapsado, verbose_=1)
        nombre_archivo = noche + "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}.fits".format(x1, x2, y1, y2,
                                                                                                 crpix1, crpix2)

        mostrarresultados(['N', 'Crpix2', 'Crpix1', 'A', 'B', '-1'],
                          [len(indice_seccion_[indice_seccion_ == seccion]), crpix2, crpix1,
                           naxis1_expected, naxis2_expected, nombre_archivo],
                          titulo='Bias Realizado')

        if cabecera is None:
            raise ValueError('No se ha creado correctamente la cabecera')

        masterbias_header = cabecera.copy()
        # if masterbias_header['BLANK']:
        #     del masterbias_header['BLANK']

        masterbias_final = fits.PrimaryHDU(master_bias_colapsado.astype(np.float), masterbias_header)

        masterbias_final.writeto(dir_bias_ + nombre_archivo, overwrite=True)

        if verbose >= 1:
            coord_lim = ImP.limites_imagen(*coordenadas_dibujo)
            ImP.imgdibujar(master_bias_colapsado, *coordenadas_dibujo, *coord_lim, verbose_=1)

        if interactive:
            # plt.show()
            input("Press Enter to continue...")


def realizar_master_biases(lista_noches, dir_listas, dir_datos, dir_bias, verbose, interactive, recortar):
    i_noche = 0
    for noche in lista_noches:
        i_noche += 1
        print('=== NOCHE ' + noche + ' - (' + str(i_noche) + '/' + str(len(lista_noches)) + ') ===')

        lista_bias = leer_lista(dir_listas + noche + '/' + 'LBias.csv')

        secciones, secciones_unicas, secciones_count, indice_seccion, bin_secciones, _ = crear_lista_unicos(
            dir_datos, noche, lista_bias, cabecera='CCDSEC', binning=True
        )

        # Variables:
        # secciones_unicas: lista de STR con las diferentes configuraciones de CCD que se usan
        # secciones_count: cuantas veces aparecen estas configuraciones
        # secciones: lista completa de STR de las secciones. se usa de apoyo. Se puede borrar despues
        # indice_seccion: INT con cual de las secciones pertenecen las calibraciones
        # coordenadas_secciones: coordenadas de las dierentes secciones
        # size-secciones: tamanyo de las imagenes en cada seccion

        coordenadas_secciones = np.zeros((len(secciones_unicas), 4), dtype=int)
        for i in range(len(secciones_unicas)):
            coordenadas_unicas = sacar_coordenadas_ccd(secciones_unicas[i])
            coordenadas_secciones[i, :] = [*coordenadas_unicas]

        # for k in range(len(secciones_unicas)):
        #     add_to_file('BiasSecciones.csv', secciones_unicas[k] + ';' + str(secciones_count[k]) + '\n')

        juntar_imagenes_bias(noche, secciones_unicas, coordenadas_secciones, secciones_count,
                             indice_seccion, bin_secciones, dir_bias, dir_datos, lista_bias,
                             verbose=verbose, interactive=interactive, recortar=recortar)


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################


def juntar_imagenes_flats(noche, secciones_unicas_, coordenadas_secciones_, indice_seccion_,
                          dir_bias_, dir_datos_, dir_flats_, lista_flats_, lista_noches, lista_bias,
                          verbose=0, interactive=False):

    # Cargamos Diccionario de Filtros
    dic_filtro = cargar_json()

    # Separamos por secciones para cada noche
    for seccion in range(len(secciones_unicas_)):
        print('seccion: ' + str(seccion))
        coordenadas_dibujo = sacar_coordenadas_2(coordenadas_secciones_, seccion)
        x1, x2, y1, y2 = deshacer_tupla_coord(coordenadas_dibujo)

        # Creamos una lista con todos aquellos que coinciden el indice
        lista_coincide = []
        for imagen in range(len(lista_flats_)):
            if indice_seccion_[imagen] == seccion:
                lista_coincide.append(lista_flats_[imagen])

        filtros, filtros_unicos, filtros_count, indice_filtro, _, nombres_filtros = crear_lista_unicos(
                                                                                      dir_datos_, noche,
                                                                                      lista_coincide,
                                                                                      cabecera='INSFLID',
                                                                                      binning=False, nombre_filtro=True
                                                                                      )

        # Creamos nombres_unicos
        nombres_unicos = []
        for i in range(len(filtros_unicos)):
            for j in range(len(filtros)):
                if filtros_unicos[i] == filtros[j]:
                    nombres_unicos.append(nombres_filtros[j])
                    break

        # Para cada elemento de nombres únicos, buscamos cual es el índice en el diccionario de filtros
        nombre_diccionario = []
        for i in nombres_unicos:
            nombre_diccionario.append(leer_diccionario(i, dic_filtro))

        # Esto sólo será necesario si NO hay grisma, luego miramos cual es el grima que le corresponde
        numero_grisma = fits.open(dir_datos_ + noche + '/' + lista_coincide[0])[0]\
            .header['insgrid']\
            .replace(' ', '0')[6:8]
        if numero_grisma == 11:  # Si numero_grisma==11, entonces no hay grisma de por medio
            free_grisma = True
        else:
            free_grisma = False

        # Queremos saber cual es el número del filtro usado (quizá se pueda quitar pues ahora se mira el nombre)
        numero_filtro = np.zeros(len(filtros_unicos), dtype=int)
        p = 0
        for i in filtros_unicos:
            numero_filtro[p] = int(i[-2:].strip())
            p += 1

        # Ahora que hemos separado en función de la sección, Separamos por los filtros
        for filtro in range(len(filtros_unicos)):
            lista_actual = []
            binning2 = []
            binning1 = []

            # Como no siempre mantienen constante los binnings, separamos por binnings
            for i in range(len(lista_coincide)):
                if indice_filtro[i] == filtro:
                    lista_actual.append(lista_coincide[i])
                    binning2.append(int(fits.open(dir_datos_ + noche + '/' + lista_coincide[i])[0].header['crpix2']))
                    binning1.append(int(fits.open(dir_datos_ + noche + '/' + lista_coincide[i])[0].header['crpix1']))

            # bin2_unique, bin2_count = np.unique(binning2, return_counts=True)
            bin1_unique, bin1_count = np.unique(binning1, return_counts=True)

            siempre_coincide_binning = True
            for i in range(len(lista_actual)):
                if not binning2 == binning1:
                    siempre_coincide_binning = False
                    break
            # solo_un_caso_de_binning = all([len(bin1_unique) == 1, len(bin2_unique) == 1])
            # if solo_un_caso_de_binning:  ############################################################################

            # Si los binnings coinciden
            master_flats_colapsado = None

            if siempre_coincide_binning:
                for binning in bin1_unique:
                    lista_actual_bin = []

                    # Ahora tenemos una lista con binnings individuales
                    for i in range(len(lista_actual)):
                        if binning1[i] == binning:
                            lista_actual_bin.append(lista_actual[i])

                    # Como tienen un solo bin y coinciden en los ejes (no se estira la imagen),
                    # no nos preocupamos y generamos las imagenes
                    crpix2 = binning
                    crpix1 = binning
                    naxis1_expected = int((coordenadas_dibujo[3] - coordenadas_dibujo[2] + 1) / crpix1)
                    naxis2_expected = int((coordenadas_dibujo[1] - coordenadas_dibujo[0] + 1) / crpix2)
                    if (coordenadas_dibujo[3] - coordenadas_dibujo[2] + 1) % crpix1 != 0:
                        naxis1_expected += 1
                    if (coordenadas_dibujo[1] - coordenadas_dibujo[0] + 1) % crpix2 != 0:
                        naxis2_expected += 1

                    master_flats = np.zeros((filtros_count[filtro], naxis1_expected, naxis2_expected), dtype=float)

                    # Generamos máscara
                    center_real = [1075, 1040]
                    center = [center_real[0] - x1, center_real[1] - y1]
                    radius = 809  # con 810 tiene un pixel de borde
                    mask = create_circular_mask(naxis1_expected, naxis2_expected, center=center, radius=radius)

                    # ToDo: Guardas las máscaras en una carpeta para que no haya que volver a hacerlas
                    # ToDo: Cronometrar crear la máscara con respecto a cargarla de disco

                    mostrarresultados(['N', 'Crpix2', 'Crpix1', 'A', 'B'],
                                      [len(lista_actual_bin), crpix2, crpix1,
                                       naxis1_expected, naxis2_expected],
                                      titulo='Filtro ' + str(filtros_count[filtro]))

                    # Leemos para cada elemento de imágenes sus valores y los vamos acoplando
                    indice0 = 0
                    cabecera = fits.open(dir_datos_ + noche + '/' + lista_actual_bin[0])[0].header
                    for imagen in range(len(lista_actual_bin)):
                        image_file = dir_datos_ + noche + '/' + lista_actual_bin[imagen]
                        image_data = fits.getdata(image_file, ext=0)

                        if image_data[:, :].shape == master_flats[indice0, :, :].shape:
                            master_flats[indice0, :, :] = image_data[:, :]                  # Juntar
                        else:
                            print('hay un problema con el tamanyo de las imagenes')
                            input('Pausado. Pulsa Enter para continuar...')
                        indice0 += 1

                    observador = cabecera['OBSERVER']
                    existe, bias_asociado_nombre, bias_asociado = obtener_bias(dir_bias_, noche, lista_noches,
                                                                               lista_bias, observador,
                                                                               x1, x2, y1, y2,
                                                                               crpix1, crpix2)

                    for i in range(master_flats.shape[0]):
                        master_flats[i, :, :] = master_flats[i, :, :] - bias_asociado

                    # Normalizamos
                    valor_medio = np.zeros(master_flats.shape[0], dtype=float)
                    for i in range(master_flats.shape[0]):

                        # Si tiene grisma, usamos una máscara para calcular la mediana
                        if free_grisma:
                            x = ma.masked_array(master_flats[i, :, :], ~mask)
                            valor_medio[i] = ma.median(x)
                        else:
                            valor_medio[i] = np.median(master_flats[i, :, :])

                        if valor_medio[i] <= 0:
                            # ToDo: Comprobar decisión tomada: Si la mediana para cada flat sale <0,
                            #  tomar el valor absoluto y si sale exactamente 0 tomar 1
                            # raise ValueError('El valor de la mediana no debe ser negativo ni 0.')

                            if valor_medio[i] == 0:
                                valor_medio[i] = 1
                            else:
                                valor_medio[i] = abs(valor_medio[i])

                        # Después ya dividimos
                        master_flats[i, :, :] = np.true_divide(master_flats[i, :, :], valor_medio[i])

                    valor_medio2 = np.zeros(master_flats.shape[0], dtype=float)
                    for i in range(master_flats.shape[0]):
                        valor_medio2[i] = np.mean(master_flats[i, :, :], dtype=float)

                    # Colapsamos
                    master_flats_colapsado = np.median(master_flats, axis=0)

                    # las esquinas las definimos como 1 para evitar problemas de dividir por 0
                    if free_grisma:
                        master_flats_colapsado[~mask] = 1

                    # plt.imshow(master_flats_colapsado)
                    # plt.show()
                    # ImP.imgdibujar(master_flats_colapsado)

                    # Generamos el nombre del fichero del flat que se va a generar
                    # si tiene numero_filtro[filtro] da la posicion del filtro
                    # si pone nombre_diccionario[filtro] da el numero del diccionario del filtro
                    nombre_archivo = noche +\
                        "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}-F{6:03d}.fits"\
                        .format(x1, x2, y1, y2, crpix1, crpix2, int(nombre_diccionario[filtro]))

                    masterflats_header = cabecera.copy()
                    if masterflats_header['BLANK']:
                        del masterflats_header['BLANK']

                    if master_flats_colapsado is None:
                        raise ValueError('No se ha creado correctamente el flat colapsado')
                    masterflats_final = fits.PrimaryHDU(master_flats_colapsado.astype(np.float), masterflats_header)

                    masterflats_final.writeto(dir_flats_ + nombre_archivo, overwrite=True)

            else:
                raise ValueError('No Coinciden los binning en ambos ejes')

            if verbose >= 1:
                coord_lim = ImP.limites_imagen(*coordenadas_dibujo)
                ImP.imgdibujar(master_flats_colapsado, *coordenadas_dibujo, *coord_lim, verbose_=1)

            if interactive:
                input("Press Enter to continue...")


def realizar_master_flats(lista_noches, lista_bias, dir_listas, dir_datos, dir_bias, dir_flats,
                          verbose, interactive):

    """
    Función maestra para generar los master flats.

    :param lista_noches:
    :param lista_bias:
    :param dir_listas:
    :param dir_datos:
    :param dir_bias:
    :param dir_flats:
    :param verbose:
    :param interactive:
    :return:
    """

    i_noche = 0
    for noche in lista_noches[:]:
        i_noche += 1
        print('=== NOCHE ' + noche + ' - (' + str(i_noche) + '/' + str(len(lista_noches)) + ') ===')

        lista_flats = leer_lista(dir_listas + noche + '/' + 'LFlat.csv')

        secciones, secciones_unicas, secciones_count, indice_seccion, bin_secciones, _ = crear_lista_unicos(
            dir_datos, noche, lista_flats, cabecera='CCDSEC', binning=True
        )

        # Variables:
        # secciones_unicas: lista de STR con las diferentes configuraciones de CCD que se usan
        # secciones_count: cuantas veces aparecen estas configuraciones
        # secciones: lista completa de STR de las secciones. se usa de apoyo. Se puede borrar despues
        # indice_seccion: INT con cual de las secciones pertenecen las calibraciones
        # coordenadas_secciones: coordenadas de las dierentes secciones
        # size-secciones: tamanyo de las imagenes en cada seccion

        coordenadas_secciones = np.zeros((len(secciones_unicas), 4), dtype=int)
        for i in range(len(secciones_unicas)):
            coordenadas_unicas = sacar_coordenadas_ccd(secciones_unicas[i])
            coordenadas_secciones[i, :] = [*coordenadas_unicas]

        juntar_imagenes_flats(noche, secciones_unicas, coordenadas_secciones, indice_seccion,
                              dir_bias, dir_datos, dir_flats, lista_flats, lista_noches, lista_bias,
                              verbose=verbose, interactive=interactive)


def realizar_reduccion(lista_noches, lista_bias, lista_flats, dir_listas, dir_datos, dir_bias, dir_flats, dir_reducc,
                       verbose, interactive):
    # Cargar listas
    dic_filtro = cargar_json()
    no_existen = []

    for noche in lista_noches:
        print(noche)
        if noche not in os.listdir(dir_reducc):
            os.mkdir(dir_reducc + noche + '/')

        lista_ciencia = leer_lista(dir_listas + noche + '/' + 'SCI.csv')
        no_existe = (0, 0, len(lista_ciencia))
        secciones, secciones_unicas, secciones_count, indice_seccion, bin_secciones, nombres_filtros = crear_lista_unicos(
            dir_datos, noche, lista_ciencia, cabecera='CCDSEC', binning=True, nombre_filtro=True
        )

        # Variables:
        # secciones_unicas: lista de STR con las diferentes configuraciones de CCD que se usan
        # secciones_count: cuantas veces aparecen estas configuraciones
        # secciones: lista completa de STR de las secciones. se usa de apoyo. Se puede borrar despues
        # indice_seccion: INT con cual de las secciones pertenecen las calibraciones
        # coordenadas_secciones: coordenadas de las dierentes secciones
        # size-secciones: tamanyo de las imagenes en cada seccion

        # print('secciones')
        # print(secciones)
        #
        # print('secciones_unicas')
        # print(secciones_unicas)
        #
        # print('secciones_count')
        # print(secciones_count)
        #
        # print('indice_seccion')
        # print(indice_seccion)
        #
        # print('bin_secciones')
        # print(bin_secciones)
        #
        # print('nombres_filtros')
        # print(nombres_filtros)

        # input('espera')

        for imagen in range(len(lista_ciencia)):
            # print('Seccion', sacar_coordenadas_ccd(secciones_unicas[indice_seccion[imagen]]))
            # print('Binning', bin_secciones[indice_seccion[imagen]])
            # print('nombre_filtro', nombres_filtros[imagen])
            # print('ID_filtro', dic_filtro[nombres_filtros[imagen]])

            # Obtenemos las coordenadas
            coordenadas = sacar_coordenadas_ccd(secciones_unicas[indice_seccion[imagen]])
            x1, x2, y1, y2 = deshacer_tupla_coord(coordenadas)

            # Optenemos el binning
            binning = bin_secciones[indice_seccion[imagen]]

            # Obtenemos el nombre del filtro y su ID_filtro
            nombrefiltro = nombres_filtros[imagen]
            id_filtro = leer_diccionario(nombrefiltro, dic_filtro)

            # Buscamos el bias
            nombre_bias = dir_bias + noche + "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}.fits"\
                .format(x1, x2, y1, y2, binning[0], binning[1])
            bias_buscado = obtener_bias(dir_bias, noche, lista_noches, lista_bias,
                                        x1, x2, y1, y2, binning[0], binning[1])

            if not bias_buscado[0]:
                print('No Existe el bias {0}'.format(nombre_bias))
                no_existe = (no_existe[0] + 1, no_existe[1], no_existe[2])

            # Buscamos el flat
            nombre_flat = dir_flats + noche + "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}-F{6:03d}.fits"\
                .format(x1, x2, y1, y2, binning[0], binning[1], int(id_filtro))
            if nombre_flat not in lista_flats:
                print('No Existe el flat {0}'.format(nombre_flat))
                no_existe = (no_existe[0], no_existe[1] + 1, no_existe[2])
            # ToDo: Hay que buscar el flat para los días anteriores y/o posteriores (igual que con los biases).
            #  Probablemente la forma más sencilla sea generando una nueva función que busque entre las listas de flats.

            # Reducimos la imagen
            # (base-bias)/flat

            # Guardamos la imagen

        # Al final de cada noche se hace el recuento
        no_existen.append(no_existe)

    mostrarresultados(lista_noches, no_existen)
    no_existen_2 = [0, 0]
    for i in range(len(no_existen)):
        no_existen_2[0] += no_existen[i][0]
        no_existen_2[1] += no_existen[i][1]

    print('Biases que no ha funcionado: ', no_existen_2[0], 'Flats que no ha funcionado: ', no_existen_2[1])



def main():

    # ---------------Valores por defecto-------------------------------------------
    default_dir_datos = '/media/enrique/TOSHIBA EXT/CAHA/CAFOS2017/'
    default_dir_bias = '/media/enrique/TOSHIBA EXT/CAHA/Biases2/'
    default_dir_listas = '/media/enrique/TOSHIBA EXT/CAHA/Listas/'
    default_dir_flats = '/media/enrique/TOSHIBA EXT/CAHA/Flats/'
    default_dir_reduccion = '/media/enrique/TOSHIBA EXT/CAHA/Reduccion/'
    desc_bias = ['bias', 'Bias', 'BIAS']
    desc_flats = ['flats', 'FLATS', 'FLAT', 'Flats', 'Flat', 'flat', 'Skyflat', 'SDkyflat']
    desc_arc = ['arc', 'ARC']
    # -----------------------------------------------------------------------------

    parser = argparse.ArgumentParser(description="Bias and Flat calibration of CAFOS images")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="store_true")
    group.add_argument("-q", "--quiet", action="store_true")
    parser.add_argument("-db", "--dir_bias", default=default_dir_bias, type=str, help='Bias Directory')
    parser.add_argument("-df", "--dir_flats", default=default_dir_flats, type=str, help='Flats Directory')
    parser.add_argument("-dd", "--dir_datos", default=default_dir_datos, type=str, help='Data Directory')
    parser.add_argument("-dl", "--dir_listas", default=default_dir_listas, type=str, help='Lists Directory')
    parser.add_argument("-de", "--dir_reducc", default=default_dir_reduccion, type=str, help='Reducction Directory')
    parser.add_argument('--cmap', type=str, help="Colormap", default='hot')
    parser.add_argument("-i", "--interactive", action="store_true")
    parser.add_argument("--recortar", action="store_true", help="Activar el recorte de imagenes")
    parser.add_argument("-nb", "--nobias", action="store_false",
                        help="No realizar los Master Bias (Por Defecto se generan).")
    parser.add_argument("-nf", "--noflat", action="store_false",
                        help="No realizar los Master Flat (Por Defecto se generan).")
    parser.add_argument("-nr", "--noreducc", action="store_false",
                        help="No realizar la reduccion (Por Defecto se generan).")
    parser.add_argument("--calysci", action="store_false",
                        help="Usar cuando los archivos no tienen '-cal-' y '-sci-'"
                             + "en el nombre para diferenciar entre calibración y ciencia.")
    args = parser.parse_args()

    # Creamos una lista de las noches disponibles
    lista_noches = os.listdir(args.dir_datos)
    tiempo_inicio = time.time()

    # Separamos entre calibración y ciencia
    crear_listas_cal_y_sci(lista_noches, args.dir_listas, args.dir_datos, desc_bias, desc_flats, desc_arc,
                           args.verbose, args.calysci)
    tiempo_listas = time.time()

    print(args.nobias, args.noflat, args.noreducc)

    # Creamos los Master Biases
    if args.nobias:
        realizar_master_biases(lista_noches, args.dir_listas, args.dir_datos, args.dir_bias,
                               args.verbose, args.interactive, args.recortar)
    lista_bias = conseguir_listas_archivos(args.dir_bias)
    tiempo_biases = time.time()

    # Creamos los Master Flats
    if args.noflat:
        realizar_master_flats(lista_noches, lista_bias, args.dir_listas, args.dir_datos, args.dir_bias, args.dir_flats,
                              args.verbose, args.interactive)
    lista_flats = conseguir_listas_archivos(args.dir_flats)
    tiempo_flats = time.time()

    # Juntamos todos los procesos y relizamos la reducción
    if args.noreducc:
        realizar_reduccion(lista_noches, lista_bias, lista_flats,
                           args.dir_listas, args.dir_datos, args.dir_bias, args.dir_flats, args.dir_reducc,
                           args.verbose, args.interactive)
    tiempo_reducc = time.time()

    # Mostramos resultados de ambos procesos
    mostrarresultados(['Tiempo Listas', 'Tiempo Master Bias', 'Tiempo Master Flats', 'Tiempo Reduccion', 'Tiempo Total',
                       'Cuantos Biases', 'Cuantos Flats'],
                      [round(tiempo_listas - tiempo_inicio, 2), round(tiempo_biases - tiempo_listas, 2),
                       round(tiempo_flats - tiempo_biases, 2), round(tiempo_reducc - tiempo_flats, 2),
                       round(tiempo_reducc - tiempo_listas, 2),
                       len(os.listdir(args.dir_bias)), len(os.listdir(args.dir_flats))],
                      titulo='Tiempo Ejecucion')


if __name__ == "__main__":

    main()
