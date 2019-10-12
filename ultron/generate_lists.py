import os
import csv
from astropy.io import fits
import numpy as np
from .Salida_limpia import stdrobust, mostrarresultados
from .auxiliary_functions import save_file_csv


def checking(file_, descriptor, descriptor2=None, descriptor3=None, verbose=False):
    if any([descriptor2, descriptor3]):
        descriptortemp = descriptor + descriptor2 + descriptor3
        coincides = True
        for texto in descriptortemp:
            if texto in fits.open(file_)[0].header['OBJECT']:
                # Some coincide without being science images
                coincides = False
                break
    else:
        coincides = False
        for texto in descriptor:
            if texto in fits.open(file_)[0].header['OBJECT']:
                # Both are the same
                coincides = True
                break
    if verbose:
        print(file_, fits.open(file_)[0].header['OBJECT'], fits.open(file_)[0].header['imagetyp'], coincides)

    # Maybe add something that checks for test

    return coincides


def create_list_cal_and_sci(lista_nights_, dir_lists_, dir_data_, desc_bias, desc_flats, desc_arc, verbose, calysci):
    i = 0
    for night in lista_nights_:
        i += 1
        if night not in os.listdir(dir_lists_):
            os.mkdir(dir_lists_ + night + '/')

            path_ = dir_data_ + night + '/'
            l_bias, l_flat, l_arc, l_ciencia, l_archivos, l_falla = file_list_2(path_,
                                                                                desc_bias,
                                                                                desc_flats,
                                                                                desc_arc,
                                                                                verbose,
                                                                                calysci)

            mostrarresultados(['Bias', 'Flat', 'Arc', 'Ciencia', 'Falla'],
                              [len(l_bias), len(l_flat), len(l_arc), len(l_ciencia), len(l_falla)],
                              titulo=night, contador=i, valor_max=len(lista_nights_))

            # save_file_csv(dir_lista_ + night + '/' + 'CAL.csv', cal)
            save_file_csv(dir_lists_ + night + '/' + 'SCI.csv', l_ciencia)
            save_file_csv(dir_lists_ + night + '/' + 'ARC.csv', l_archivos)
            save_file_csv(dir_lists_ + night + '/' + 'LBias.csv', l_bias)
            save_file_csv(dir_lists_ + night + '/' + 'LFlat.csv', l_flat)
            save_file_csv(dir_lists_ + night + '/' + 'LArc.csv', l_arc)
            save_file_csv(dir_lists_ + night + '/' + 'Errores.csv', l_falla)


def create_unique_list(dir_data, night, stuff_list, header, binning=False, filter_name=False):

    """
    Given a list of parameter, generates a series of lists:
        - List of files
        - List of unique values
        - For each value, how many times it appears
        - For every element of the list, generates a value of which 'unique value' appears
        - How many unique values does it have
        - Name of filter, if applicable

    :param dir_data:
    :param night:
    :param stuff_list:
    :param header:
    :param binning:
    :param filter_name:
    :return:
    """
    lista = []

    for image in stuff_list:
        lista.append(fits.open(dir_data + night + '/' + image)[0].header[header])
    unique_list, count_list = np.unique(lista, return_counts=True)  # Counting each list

    bin_sections = -1 * np.ones((len(unique_list), 2), dtype=int)
    indx_stuff = np.zeros(len(stuff_list), dtype=int)
    name_filter = []

    for i in range(len(stuff_list)):
        for j in range(len(unique_list)):
            if lista[i] == unique_list[j]:
                indx_stuff[i] = j
                if binning:
                    bin_sections[j, 0] = int(fits.open(dir_data + night + '/' + stuff_list[i])[0].header['ccdbinX'])
                    bin_sections[j, 1] = int(fits.open(dir_data + night + '/' + stuff_list[i])[0].header['ccdbinY'])
                if filter_name:
                    name_filter.append(fits.open(dir_data + night + '/' + stuff_list[i])[0].header['INSFLNAM'])
                break

    return lista, unique_list, count_list, indx_stuff, bin_sections, name_filter


def file_list_2(path_, desc_bias, desc_flats, desc_arc, verbose=False, calysci=True):
    lista_cal = []
    lista_sci = []
    lista_misc = []
    lista_bias = []
    lista_flat = []
    lista_science = []
    lista_arc = []
    lista_else = []
    lista_files = []
    lista_wrong = []

    # Listing every file in the folder
    for file in os.listdir(path_):
        if file.endswith(".fits"):
            if calysci:
                if '-cal-' in file:
                    lista_cal.append(file)
                elif '-sci-' in file:
                    lista_sci.append(file)
                else:
                    lista_misc.append(file)
            lista_files.append(os.path.join(path_, file))

            # Splitting by IMAGETYP
            types = fits.open(path_ + file)[0].header['IMAGETYP'].strip()
            if not fits.open(path_ + file)[0].header['OBJECT'] == 'Test':  # Check it's not a test file
                if types == 'BIAS' or types == 'bias':
                    coincides = checking(path_ + file, desc_bias, verbose=verbose)
                    if coincides:
                        lista_bias.append(file)
                    else:
                        lista_wrong.append(file)

                elif types == 'flat':
                    coincides = checking(path_ + file, desc_flats, verbose=verbose)
                    if coincides:
                        lista_flat.append(file)
                    else:
                        lista_wrong.append(file)

                elif types == 'arc':
                    coincides = checking(path_ + file, desc_arc, verbose=verbose)
                    if coincides:
                        lista_arc.append(file)
                    else:
                        lista_wrong.append(file)

                elif types == 'science':
                    coincides = checking(path_ + file, desc_bias, desc_flats, desc_arc, verbose=verbose)
                    if coincides:
                        lista_science.append(file)
                    else:
                        lista_wrong.append(file)

                else:
                    lista_else.append(file)

    for file in lista_wrong:

        if calysci:
            if file in lista_sci:
                lista_science.append(file)

        probable = most_probable_image(path_ + file)

        if probable == 0:
            lista_bias.append(file)
        elif probable == 1:
            lista_arc.append(file)
        elif probable == 2:
            lista_flat.append(file)

    return lista_bias, lista_flat, lista_arc, lista_science, lista_files, lista_wrong


def most_probable_image(archivo):
    image_data = fits.getdata(archivo, ext=0)
    v_median = np.median(image_data)
    v_sd = stdrobust(image_data)
    if v_sd < 50 and v_median < 900:
        probable = 0  # Bias
    elif v_sd < 500:
        probable = 1  # Arc
    else:
        probable = 2  # Flat

    return probable


def obtain_files_lists(path_):
    file_list = []
    for file in os.listdir(path_):
        if file.endswith(".fits"):
            file_list.append(os.path.join(path_, file))
    return file_list


def read_list(archivo):
    with open(archivo, 'rt') as f:
        reader = csv.reader(f, delimiter=',')
        your_list = list(reader)
    your_list = [item for sublist in your_list for item in sublist]
    return your_list
