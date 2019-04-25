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


def stdrobusta(array):
    return 0.7413 * (np.quantile(array.reshape(-1,), 0.75) - np.quantile(array.reshape(-1,), 0.25))


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
    df = pd.DataFrame(columns=['Medio', 'Mediana', 'Relacion', 'Std', 'Stdrob', 'Temp', 'Fecha'])

    i = 0
    for imagen in lista_noches:
        hdul = fits.open(args.dir_bias + imagen)
        image_header = hdul[0].header
        image_data = fits.getdata(args.dir_bias + imagen, ext=0)
        hdul.close()
        valor_std = np.std(image_data)
        valor_stdrob = stdrobusta(image_data)
        valor_medio = np.mean(image_data)
        valor_mediana = np.median(image_data)
        rel_medio_mediana = (valor_medio-valor_mediana)/valor_mediana*100
        temperatura = image_header['CCDTEMP']
        # fechajul = image_header['MJD-OBS']
        fechanorm = image_header['DATE']
        your_date = datetime.datetime.strptime(fechanorm, "20%y-%m-%dT%H:%M:%S")
        i += 1
        lista_imagen = [valor_medio, valor_mediana, rel_medio_mediana, valor_std, valor_stdrob, temperatura, your_date]
        df.loc[i] = np.asarray(lista_imagen)

        # print(valor_medio, valor_varia, temperatura)
        # input("continuar")
    df = df.sort_values(by='Fecha', ascending=True).reset_index()
    df = df.drop(df.index[13])

    print(df)
    eje_x = df.Fecha.values
    medianas = df.Medio.values
    stdrob = df.Stdrob.values
    temp = df.Temp.values

    #######################################################
    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('fecha')
    ax1.set_ylabel('Cuentas', color=color)
    # ax1.scatter(eje_x, medianas, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.errorbar(eje_x, medianas, stdrob, color=color, fmt='o')

    # ax1.plot(eje_x, medias + stdrob, color=color)
    # ax1.plot(eje_x, medias - stdrob, color=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:green'
    ax2.set_ylabel('Temperatura', color=color)  # we already handled the x-label with ax1
    ax2.scatter(eje_x, temp, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.title('Evolucion del Bias a lo largo del año y la temperatura')
    plt.show()






if __name__ == "__main__":

    main()


# matplotlilb plot_date
# Date Demo Rrule
