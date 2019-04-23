from astropy.io import fits
from Salida_limpia import mostrarresultados
import numpy as np
import IMGPlot as ImP
import pandas as pd
import matplotlib.pyplot as plt
import os
import csv
from numpy import genfromtxt
import argparse
import time
import warnings
import datetime
# your_date = datetime.datetime.strptime("31/12/2016", "%d/%m/%Y")


def transformar_a_int(fecha):
    ye = fecha.year
    mo = fecha.month
    da = fecha.day
    ho = fecha.hour
    mi = fecha.minute
    se = fecha.second

    numero = ho/24 + mi/(24*60) + se/(24*3600) + 365*ye + mo


def main():
    # --------------- Valores por defecto -----------------------------------------
    default_dir_datos = 'CAFOS2017/'
    default_dir_bias = '/media/enrique/TOSHIBA EXT/CAHA/Biases2/'
    default_dir_listas = 'Listas/'
    # -----------------------------------------------------------------------------

    parser = argparse.ArgumentParser(description="Bias and Flat calibration of CAFOS images")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="store_true")
    group.add_argument("-q", "--quiet", action="store_true")
    parser.add_argument("-db", "--dir_bias", default=default_dir_bias, type=str, help='Bias Directory')
    parser.add_argument("-dd", "--dir_datos", default=default_dir_datos, type=str, help='Data Directory')
    parser.add_argument("-dl", "--dir_listas", default=default_dir_listas, type=str, help='Lists Directory')
    parser.add_argument('--cmap', type=str, help="Colormap", default='hot')
    parser.add_argument("-i", "--interactive", action="store_true")
    parser.add_argument("--recortar", action="store_true")
    args = parser.parse_args()

    lista_noches = os.listdir(args.dir_bias)
    df = pd.DataFrame(columns=['FechaJ', 'Medio', 'Var', 'Temp', 'Fecha'])

    dfq = pd.DataFrame(columns=['lib', 'qty1', 'qty2'])
    for i in range(5):
        dfq.loc[i] = [np.random.randint(-1, 1) for n in range(3)]



    lista_fecha = []
    i = 0
    for imagen in lista_noches:
        # print(imagen)
        hdul = fits.open(args.dir_bias + imagen)
        image_header = hdul[0].header
        image_data = fits.getdata(args.dir_bias + imagen, ext=0)
        hdul.close()
        plt.imshow(image_data)
        plt.show()
        print(image_data)
        valor_varia = np.var(image_data)
        valor_medio = np.mean(image_data)
        temperatura = image_header['CCDTEMP']
        fechajul = image_header['MJD-OBS']
        fechanorm = image_header['DATE']
        lista_imagen = [fechajul, valor_medio, valor_varia, temperatura, fechanorm]
        df.loc[i] = np.asarray(lista_imagen)
        i += 1
        your_date = datetime.datetime.strptime(fechanorm, "20%y-%m-%dT%H:%M:%S")
        lista_fecha.append(your_date)
        print(your_date.year, your_date.day, your_date.month, your_date.hour, your_date.minute, your_date.second)

        # print(valor_medio, valor_varia, temperatura)
        # input("continuar")

    eje_x = df.Fecha.values
    eje_y = df.Medio.values
    print(eje_y)

    plt.plot(lista_fecha, eje_y)
    plt.show()


if __name__ == "__main__":

    main()


# matplotlilb plot_date
# Date Demo Rrule
