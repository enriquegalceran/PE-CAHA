# -*- coding: utf-8 -*-

"""
    Proyecto de Unidad de Limpieza y Tratamiento de Resultados Observacionales Nativo (Proyecto U.L.T.R.O.N.)

    Recoge los resultados de observaciones del instrumento CAFOS y los reduce correctamente.
    Desarrollado por Enrique Galceran García
"""

import pandas as pd
import os
import argparse
import time

from .Salida_limpia import mostrarresultados
from .generate_lists import obtain_files_lists, create_list_cal_and_sci
from .bias import make_master_bias
from .flat import make_master_flat
from .reduction import reducing_images, decidir_repetir_calculos


def main():

    # ---------------Valores por defecto-------------------------------------------
    default_dir_datos = '/media/enrique/TOSHIBA EXT/DataCAHA/CAFOS2017/'
    default_dir_bias = '/media/enrique/TOSHIBA EXT/CAHA/Biases/'
    default_dir_listas = '/media/enrique/TOSHIBA EXT/CAHA/Listas/'
    default_dir_flats = '/media/enrique/TOSHIBA EXT/CAHA/Flats/'
    default_dir_reduccion = '/media/enrique/TOSHIBA EXT/CAHA/Reduccion/'
    default_dir_dataframe = ''
    desc_bias = ['bias', 'Bias', 'BIAS']
    desc_flats = ['flats', 'FLATS', 'FLAT', 'Flats', 'Flat', 'flat', 'Skyflat', 'SDkyflat']
    desc_arc = ['arc', 'ARC']
    # -----------------------------------------------------------------------------

    parser = argparse.ArgumentParser(description="Bias and Flat calibration of CAFOS images")
    group = parser.add_mutually_exclusive_group()
    grupo2 = parser.add_mutually_exclusive_group()
    grupo3 = parser.add_mutually_exclusive_group()

    group.add_argument("-v", "--verbose", action="store_true")
    group.add_argument("-q", "--quiet", action="store_true")
    parser.add_argument("-db", "--dir_bias", default=default_dir_bias, type=str, help='Bias Directory')
    parser.add_argument("-df", "--dir_flats", default=default_dir_flats, type=str, help='Flats Directory')
    parser.add_argument("-dd", "--dir_datos", default=default_dir_datos, type=str, help='Data Directory')
    parser.add_argument("-dl", "--dir_listas", default=default_dir_listas, type=str, help='Lists Directory')
    parser.add_argument("-de", "--dir_reducc", default=default_dir_reduccion, type=str, help='Reducction Directory')
    parser.add_argument("-ddf", "--dir_dataf", default=default_dir_dataframe, type=str, help='DataFrame Directory')
    parser.add_argument("-vi", "--verboseimage", action="store_true", help="Mostrar Imagenes")
    parser.add_argument('--cmap', type=str, help="Colormap", default='hot')
    parser.add_argument("-i", "--interactive", action="store_true")
    parser.add_argument("--recortar", action="store_true", help="Activar el recorte de imagenes")
    parser.add_argument("--calysci", action="store_false",
                        help="Usar cuando los archivos no tienen '-cal-' y '-sci-'"
                             + "en el nombre para diferenciar entre calibración y ciencia.")
    parser.add_argument("-nr", "--noreducc", action="store_false", help="No realizar la reduccion.")

    grupo2.add_argument("-nb", "--nobias", action="store_true",
                        help="No realizar los Master Bias.")
    grupo2.add_argument("-sb", "--sibias", action="store_true",
                        help="Fuerza realizar los Master Bias.")

    grupo3.add_argument("-nf", "--noflat", action="store_true",
                        help="No realizar los Master Flat.")
    grupo3.add_argument("-sf", "--siflat", action="store_true",
                        help="Fuerza realizar los Master Flat.")

    args = parser.parse_args()

    if args.verbose:
        verbosidad = 2
    elif args.quiet:
        verbosidad = 0
    else:
        verbosidad = 1
    print('verbosidad', verbosidad)

    # Comprobamos si queremos/hace falta calcular los bias/flats
    print('bias:', args.nobias, args.sibias)
    realizarbias = decidir_repetir_calculos(args.nobias, args.sibias, 'bias', args.dir_dataf, args.dir_bias)
    print('flat:', args.noflat, args.siflat)
    realizarflat = decidir_repetir_calculos(args.noflat, args.siflat, 'flat', args.dir_dataf, args.dir_flats)

    # Creamos una lista de las noches disponibles
    lista_noches = os.listdir(args.dir_datos)
    lista_noches.sort()
    tiempo_inicio = time.time()

    # Separamos entre calibración y ciencia
    create_list_cal_and_sci(lista_noches, args.dir_listas, args.dir_datos, desc_bias, desc_flats, desc_arc,
                            verbosidad, args.calysci)
    tiempo_listas = time.time()

    print(args.nobias, args.noflat, args.noreducc)

    # importado_b = pd.read_csv('df_bias.csv')
    # importado_f = pd.read_csv('df_flat.csv')
    # pruebas_pandas(importado_f)

    # Creamos los Master Biases
    if realizarbias:
        df_bias = make_master_bias(lista_noches, args.dir_listas, args.dir_datos, args.dir_bias,
                                   args.interactive, args.recortar, verbose_imagen=args.verboseimage)
        numero_bias = len(os.listdir(args.dir_bias))
        print(df_bias)
        _ = df_bias.to_csv('df_bias.csv', index=None, header=True)
    else:
        df_bias = pd.read_csv('df_bias.csv')
        print('Se han importado los bias')
        numero_bias = '-'

    lista_bias = obtain_files_lists(args.dir_bias)
    tiempo_biases = time.time()

    # Creamos los Master Flats
    if realizarflat:
        df_flat = make_master_flat(lista_noches, lista_bias,
                                   args.dir_listas, args.dir_datos, args.dir_bias, args.dir_flats,
                                   verbosidad, args.interactive, verbose_imagen=args.verboseimage)
        numero_flats = len(os.listdir(args.dir_flats))
        print(df_flat)
        _ = df_flat.to_csv('df_flat.csv', index=None, header=True)
    else:
        df_flat = pd.read_csv('df_flat.csv')
        print('Se han importado los flats')
        numero_flats = '-'

    tiempo_flats = time.time()

    # Juntamos todos los procesos y relizamos la reducción
    if args.noreducc:
        numeros_reducidos = reducing_images(lista_noches, args.dir_listas, args.dir_datos, args.dir_bias,
                                            args.dir_flats, args.dir_reducc, df_bias, df_flat,
                                            verbosidad, verbose_imagen=args.verboseimage)
    else:
        numeros_reducidos = '-'
    tiempo_reducc = time.time()

    # Mostramos resultados de ambos procesos
    mostrarresultados(['Tiempo Listas', 'Tiempo Master Bias', 'Tiempo Master Flats', 'Tiempo Reduccion', 'Tiempo Total',
                       'Cuantos Biases', 'Cuantos Flats', 'Cuantos Reducidos'],
                      [round(tiempo_listas - tiempo_inicio, 2), round(tiempo_biases - tiempo_listas, 2),
                       round(tiempo_flats - tiempo_biases, 2), round(tiempo_reducc - tiempo_flats, 2),
                       round(tiempo_reducc - tiempo_listas, 2),
                       numero_bias, numero_flats, numeros_reducidos],
                      titulo='Tiempo Ejecucion')


if __name__ == "__main__":

    main()
