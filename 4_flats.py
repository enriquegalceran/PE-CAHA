from astropy.io import fits
from Salida_limpia import mostrarresultados  # , stdrobusta
import numpy as np
import numpy.ma as ma
import IMGPlot as ImP
import os
import csv
import argparse
import time
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
    json_f = json.dumps(variable)
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
    for x in range(1, longitud):
        if str_ccdsec[x] == ',':
            coma = x
            break
    for x in range(coma + 1, longitud):
        if str_ccdsec[x] == ':':
            puntos = x
            break
    for x in range(puntos + 1, longitud):
        if str_ccdsec[x] == ',':
            coma2 = x
            break

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
    return your_list


def sacar_coordenadas_2(lista, idx):
    x1 = lista[idx, 0]
    x2 = lista[idx, 1]
    y1 = lista[idx, 2]
    y2 = lista[idx, 3]
    tupla_salida = (x1, x2, y1, y2)
    return tupla_salida


def obtener_bias(dir_bias_, noche, lista_noches, lista_bias, observador, x1, x2, y1, y2, b1, b2, busqueda_max=10):

    """
    Busca en la lista de bias disponibles si exite un archivo de bias para la misma noche en la que se hizo la foto.
    Si exite, usará esa directamente. Si no existe, busca una imagen de bias diferente del mismo autor de noches
    cercanas. Si no hay, genera un bias plano.

    :param dir_bias_:
    :param noche:
    :param lista_noches:
    :param lista_bias:
    :param observador:
    :param x1:
    :param x2:
    :param y1:
    :param y2:
    :param b1:
    :param b2:
    :param busqueda_max:
    :return:
    """

    bias_asociado_nombre = dir_bias_ + noche + "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}.fits".format(x1, x2,
                                                                                                               y1, y2,
                                                                                                               b1, b2)
    otro_observador = False
    hace_falta_cambio = True
    indice_bias = None
    if bias_asociado_nombre in lista_bias:
        hace_falta_cambio = False
        bias_asociado = fits.getdata(bias_asociado_nombre, ext=0)
    else:
        exitos = []
        observadores = []
        posicion = lista_noches.index(noche)
        for i in range(1, busqueda_max):
            for mult in [-1, 1]:
                indice = i * mult
                noche = lista_noches[posicion + indice]
                bias_asociado_nuevo = dir_bias_ + noche + \
                    "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}.fits"\
                    .format(x1, x2, y1, y2, b1, b2)

                if bias_asociado_nuevo in lista_bias:
                    exitos.append(indice)
                    observadores.append(fits.open(bias_asociado_nuevo)[0].header['OBSERVER'])

        # Ahora tenemos una lista con posibles candidatos
        if all([exitos, observadores]):
            for j in range(len(observadores)):
                if observadores[j] == observador:
                    indice_bias = exitos[j]
                    break

            # Si no ha encontrado nada, cogemos el primero
            if not indice_bias:
                indice_bias = exitos[0]
                otro_observador = True
        else:
            print('No hay biases cerca')
            input('Enter para continuar...')
        if indice_bias:

            noche = lista_noches[posicion + indice_bias]
            bias_asociado_n = dir_bias_ + noche +\
                "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}.fits".format(x1, x2, y1, y2, b1, b2)
            bias_asociado = fits.getdata(bias_asociado_n, ext=0)

    if hace_falta_cambio:
        print(exitos)
        print(observadores)
    else:
        print('existe y sin problemas')

    return bias_asociado, noche, otro_observador


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

            bin2_unique, bin2_count = np.unique(binning2, return_counts=True)
            bin1_unique, bin1_count = np.unique(binning1, return_counts=True)

            siempre_coincide_binning = True
            for i in range(len(lista_actual)):
                if not binning2 == binning1:
                    siempre_coincide_binning = False
                    break
            # solo_un_caso_de_binning = all([len(bin1_unique) == 1, len(bin2_unique) == 1])
            # if solo_un_caso_de_binning:  ############################################################################

            # Si los binnings coinciden
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

                    if free_grisma:
                        center_real = [1075, 1040]
                        center = [center_real[0] - x1, center_real[1] - y1]
                        radius = 809  # con 810 tiene un pixel de borde
                        mask = create_circular_mask(naxis1_expected, naxis2_expected, center=center, radius=radius)

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
                    bias_asociado, noche_usada, otro_observador = obtener_bias(dir_bias_, noche, lista_noches,
                                                                               lista_bias, observador,
                                                                               x1, x2, y1, y2,
                                                                               crpix1, crpix2)
                    if noche_usada != noche:
                        print('Se ha usado una noche diferente: ', noche_usada)
                    if otro_observador:
                        print('La imagen es de otro observador')

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

                        # Después ya dividimos
                        master_flats[i, :, :] = np.true_divide(master_flats[i, :, :], valor_medio[i])

                    valor_medio2 = np.zeros(master_flats.shape[0], dtype=float)
                    for i in range(master_flats.shape[0]):
                        valor_medio2[i] = np.mean(master_flats[i, :, :], dtype=float)

                    # Colapsamos
                    master_flats_colapsado = np.median(master_flats, axis=0)
                    # plt.imshow(master_flats_colapsado)
                    # plt.show()
                    # ImP.imgdibujar(master_flats_colapsado)

                    # Generamos el nombre del fichero del flat que se va a generar
                    # si tiene numero_filtro[filtro] da la posicion del filtro
                    # si pone nombre_diccionario[filtro] da el numero del diccionario del filtro
                    nombre_archivo = noche +\
                        "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}-F{6:03d}.fits"\
                        .format(x1, x2, y1, y2, crpix1, crpix2, nombre_diccionario[filtro])

                    masterflats_header = cabecera.copy()
                    if masterflats_header['BLANK']:
                        del masterflats_header['BLANK']

                    masterflats_final = fits.PrimaryHDU(master_flats_colapsado.astype(np.float), masterflats_header)

                    masterflats_final.writeto(dir_flats_ + nombre_archivo, overwrite=True)

            else:
                print('no coinciden los binnings')
                input('Pausa')

            if verbose >= 1:
                coord_lim = ImP.limites_imagen(*coordenadas_dibujo)
                ImP.imgdibujar(master_flats_colapsado, *coordenadas_dibujo, *coord_lim, verbose_=1)

            if interactive:
                input("Press Enter to continue...")


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
        lista_flats = [item for sublist in lista_flats for item in sublist]  # Limpiamos la lista para poder usarla

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


def main():
    # ---------------Valores por defecto-------------------------------------------
    default_dir_datos = 'CAFOS2017/'
    default_dir_bias = 'Biases/'
    default_dir_listas = 'Listas/'
    default_dir_flats = 'Flats/'
    # -----------------------------------------------------------------------------

    parser = argparse.ArgumentParser(description="Bias and Flat calibration of CAFOS images")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="store_true")
    group.add_argument("-q", "--quiet", action="store_true")
    parser.add_argument("-db", "--dir_bias", default=default_dir_bias, type=str, help='Bias Directory')
    parser.add_argument("-df", "--dir_flats", default=default_dir_flats, type=str, help='Flats Directory')
    parser.add_argument("-dd", "--dir_datos", default=default_dir_datos, type=str, help='Data Directory')
    parser.add_argument("-dl", "--dir_listas", default=default_dir_listas, type=str, help='Lists Directory')
    parser.add_argument('--cmap', type=str, help="Colormap", default='hot')
    parser.add_argument("-i", "--interactive", action="store_true")
    args = parser.parse_args()

    lista_noches = os.listdir(args.dir_datos)
    tiempo_inicio_listas = time.time()

    lista_bias = conseguir_listas_archivos(args.dir_bias)
    realizar_master_flats(lista_noches, lista_bias, args.dir_listas, args.dir_datos, args.dir_bias, args.dir_flats,
                          args.verbose, args.interactive)
    tiempo_final = time.time()

    mostrarresultados(['Tiempo Master Bias', 'Cuantos Flats'],
                      [round(tiempo_final - tiempo_inicio_listas, 2),
                       len(os.listdir(args.dir_flats))],
                      titulo='Tiempo Ejecucion')


if __name__ == "__main__":

    main()

################################################################
# Dibujar uno de ellos

# print('DIBUJAR')
# image_file = mypath + lista_cal[0]

# hdul = fits.open(image_file)
# hdr = hdul[0].header

# #print(hdr[0:5])
# print(hdr['CCDSEC'])
# # for imagen in lista_cal:
# #     print(fits.open(mypath + imagen)[0].header['OBJECT'])


# image_data = fits.getdata(image_file, ext=0)
# plt.imshow(image_data)
# plt.show()


# Fijar nuevas variables en el header:
# hdr = hdul[0].header
# >>> hdr['targname'] = ('NGC121-a', 'the observation target')
# >>> hdr['targname']
# 'NGC121-a'
# >>> hdr.comments['targname']
# 'the observation target'

# >>> hdr.set('observer', 'Edwin Hubble')

# History y Comments anyade un valor nuevo
# >>> hdr['history'] = 'I updated this file 2/26/09'
# >>> hdr['comment'] = 'Edwin Hubble really knew his stuff'
# >>> hdr['comment'] = 'I like using HST observations'
# >>> hdr['history']
# I updated this file 2/26/09
# >>> hdr['comment']
# Edwin Hubble really knew his stuff
# I like using HST observations

# Cambiar el valor es muy sencillo, sigue como numpy
# >>> data[30:40, 10:20] = data[1, 4] = 999
# >>> photflam = hdul[1].header['photflam']
# >>> exptime = hdr['exptime']
# >>> data = data * photflam / exptime
# >>> hdul.close()


###################################################################
# Dibujar una imagen
# #image_file = get_pkg_data_filename('tutorials/FITS-images/HorseHead.fits')
# image_file = 'caf-20170224-20:20:04-cal-krek.fits'


# fits.info(image_file)

# hdu = fits.PrimaryHDU()
# print('%&%&%&%&%&%&%&%&%')
# #print(hdu.header)

# hdul = fits.open(image_file)
# hdr = hdul[0].header
# #print(hdr[34])

# #print(hdr)
# print(hdr['OBJECT'])
# print(hdr['BIASSEC'])
# print(hdr['DATASEC'])
# print(hdr['CCDSEC'])


# image_data = fits.getdata(image_file, ext=0)

# print(image_data.shape)

# plt.figure()
# plt.imshow(image_data, cmap='gray')
# plt.colorbar()
# plt.show()


# '#!/bin/bash'
# 'python 3_inicio.py nombre archivo\
# '--parametro1 valor1 \
# '--parametro2 valor2 \
# chmod +x nombre_script.sh

# https://pyformat.info
# "str{0:04d}".format(210)


# usar numpy.ma.median para la mascara
# Crear una mascara que calcule el radio con mask=x*x+y*y<=r*r
# Dividir por median
#
# Anadir B{0:02d}-{1:02d}.format(binning1,binning2)
# Hay que comprobar a ver si hay imagenes de ciencia con multiples binning tambien o no
#
# Crear un diccionario de filtros que se van a usar para el nombre. Siempre se usaran los mimso numeros para el mismo
# filtro aunque sean diferentes noches o imagenes. Ese diccionario se ampliara a lo largo del tiempo
