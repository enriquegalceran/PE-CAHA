from astropy.io import fits
from astropy.time import Time
import datetime
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd

from .auxiliary_functions import tuple2coordinates, obtain_naxis
from .bias import obtain_bias
from .coordinates import obtain_coordinates_2, obtain_coordinates_ccd
from .dictionary import read_dictionary
from .generate_lists import read_list, create_unique_list
from .IMGPlot import imgdibujar, limites_imagen
from .json_functions import load_json
from .mask_generation import create_circular_mask
from .Salida_limpia import mostrarresultados


def join_flat_images(nights, unique_sections_, sections_coordinates_, indx_section_,
                     dir_bias_, dir_data_, dir_flats_, list_flats_, list_nights_, lista_bias,
                     verbose=0, interactive=False, verbose_images=False):
    """
    Combine multiple flats to create a master flat.

    :param nights:
    :param unique_sections_:
    :param sections_coordinates_:
    :param indx_section_:
    :param dir_bias_:
    :param dir_data_:
    :param dir_flats_:
    :param list_flats_:
    :param list_nights_:
    :param lista_bias:
    :param verbose:
    :param interactive:
    :param verbose_images:
    :return:
    """

    # Load Dictionary File
    dic_filter = load_json()
    element_list = None

    # Splitting by night
    for section in range(len(unique_sections_)):
        print('section: ' + str(section))
        coordinates_image = obtain_coordinates_2(sections_coordinates_, section)
        x1, x2, y1, y2 = tuple2coordinates(coordinates_image)

        # Create a list with every one that fits the index
        list_coincide = []
        for image in range(len(list_flats_)):
            if indx_section_[image] == section:
                list_coincide.append(list_flats_[image])

        filters, unique_filters, filter_counts, indx_filter, _, name_filter = create_unique_list(
            dir_data_, nights, list_coincide, header='INSFLID', binning=False, filter_name=True
        )

        # crete unique_names
        unique_names = []
        for i in range(len(unique_filters)):
            for j in range(len(filters)):
                if unique_filters[i] == filters[j]:
                    unique_names.append(name_filter[j])
                    break

        # For each element of unique_names, we search it's intex in the dictionary
        name_dictionary = []
        for i in unique_names:
            name_dictionary.append(read_dictionary(i, dic_filter))

        # Only necessary if the ISN'T a grisma, so we search what it's grisma is being used
        grisma_number = int(fits.open(dir_data_ + nights + '/' + list_coincide[0])[0]
                            .header['insgrid'].replace(' ', '0')[6:8])

        if grisma_number == 11:  # If grisma_number==11, there is no grisma being used
            free_grisma = True
        else:
            free_grisma = False

        # We want to know what is the number of the filter used
        filter_number = np.zeros(len(unique_filters), dtype=int)
        p = 0
        for i in unique_filters:
            filter_number[p] = int(i[-2:].strip())
            p += 1

        # Separated by sections, so we separate by filter
        for filtro in range(len(unique_filters)):
            lista_actual = []
            binning2 = []
            binning1 = []

            # Binning is constant, so we split by different binning
            for i in range(len(list_coincide)):
                if indx_filter[i] == filtro:
                    lista_actual.append(list_coincide[i])
                    binning2.append(int(fits.open(dir_data_ + nights + '/' + list_coincide[i])[0].header['ccdbinX']))
                    binning1.append(int(fits.open(dir_data_ + nights + '/' + list_coincide[i])[0].header['ccdbinY']))

            bin1_unique, bin1_count = np.unique(binning1, return_counts=True)

            binning_always_coincide = True
            for i in range(len(lista_actual)):
                if not binning2 == binning1:
                    binning_always_coincide = False
                    break

            # If binning coincide
            master_flats_colapsed = None

            if binning_always_coincide:
                for binning in bin1_unique:
                    lista_actual_bin = []

                    # We have a list with individual binnings
                    for i in range(len(lista_actual)):
                        if binning1[i] == binning:
                            lista_actual_bin.append(lista_actual[i])

                    # Only has one bin and it matches, it's axis (image isn't stretched),
                    # we don't worry about it and generate the image
                    ccdbiny = binning
                    ccdbinx = binning

                    naxis1_expected, naxis2_expected = obtain_naxis(coordinates_image, ccdbinx, ccdbiny)

                    master_flats = np.zeros((filter_counts[filtro], naxis1_expected, naxis2_expected), dtype=float)

                    # Generate mask
                    real_center = [1075, 1040]
                    center = [real_center[0] - x1, real_center[1] - y1]
                    radius = 809  # con 810 tiene un pixel de borde
                    mask = create_circular_mask(naxis1_expected, naxis2_expected, center=center, radius=radius)

                    mostrarresultados(['N', 'ccdbiny', 'ccdbinx', 'A', 'B'],
                                      [len(lista_actual_bin), ccdbiny, ccdbinx,
                                       naxis1_expected, naxis2_expected],
                                      titulo='Filtro ' + str(name_dictionary[filtro]))

                    # Read for each element their values and we concatenate
                    index0 = 0
                    header_ = fits.open(dir_data_ + nights + '/' + lista_actual_bin[0])[0].header
                    for image in range(len(lista_actual_bin)):
                        image_file = dir_data_ + nights + '/' + lista_actual_bin[image]
                        image_data = fits.getdata(image_file, ext=0)

                        if image_data[:, :].shape == master_flats[index0, :, :].shape:
                            master_flats[index0, :, :] = image_data[:, :]  # Juntar
                        else:
                            print('There is a problem with the image size')
                            input('Pause. Press Enter to to continue...')
                        index0 += 1

                    exist, searched_bias_name, searched_bias = obtain_bias(
                        dir_bias_, nights, list_nights_, lista_bias, x1, x2, y1, y2, ccdbinx, ccdbiny
                    )

                    for i in range(master_flats.shape[0]):
                        master_flats[i, :, :] = master_flats[i, :, :] - searched_bias

                    # Normalize
                    median_value = np.zeros(master_flats.shape[0], dtype=float)
                    for i in range(master_flats.shape[0]):

                        # If it has a grisma, we need to use the mask
                        if free_grisma:
                            x = ma.masked_array(master_flats[i, :, :], ~mask)
                            median_value[i] = ma.median(x)
                        else:
                            median_value[i] = np.median(master_flats[i, :, :])

                        if median_value[i] <= 0:
                            if median_value[i] == 0:
                                median_value[i] = 1
                            else:
                                median_value[i] = abs(median_value[i])

                        # Divide
                        master_flats[i, :, :] = np.true_divide(master_flats[i, :, :], median_value[i])

                    mean_value2 = np.zeros(master_flats.shape[0], dtype=float)
                    for i in range(master_flats.shape[0]):
                        mean_value2[i] = np.mean(master_flats[i, :, :], dtype=float)

                    # Collapse
                    master_flats_colapsed = np.median(master_flats, axis=0)

                    # Set corners as 1 to avoid dividing by 0
                    if free_grisma:
                        master_flats_colapsed[~mask] = 1

                    if verbose_images:
                        imgdibujar(master_flats_colapsed)
                        plt.show()

                    # Generate the name of the future flat
                    flat_file_name = nights + \
                        "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}-F{6:03d}.fits" \
                        .format(x1, x2, y1, y2, ccdbinx, ccdbiny, int(name_dictionary[filtro]))

                    masterflats_header = header_.copy()
                    if masterflats_header['BLANK']:
                        del masterflats_header['BLANK']

                    now_time = datetime.datetime.now()
                    now_datetime = Time(now_time, format='datetime', scale='utc')
                    masterflats_header.add_history('Realizado Flat a partir de ' +
                                                   str(len(indx_section_[indx_section_ == section])) +
                                                   ' imagenes. | ' + str(now_datetime)[:19])
                    flat_number = 0
                    for flat_raw in range(len(list_flats_)):
                        if indx_section_[flat_raw] == section:
                            flat_number += 1
                            masterflats_header.add_history(str(flat_number) + ': ' + list_flats_[flat_raw])
                            print(str(flat_number) + ': ' + list_flats_[flat_raw])

                    if master_flats_colapsed is None:
                        raise ValueError('Collapsed flat created wrong.')
                    masterflats_final = fits.PrimaryHDU(master_flats_colapsed.astype(np.float32), masterflats_header)

                    masterflats_final.writeto(dir_flats_ + flat_file_name, overwrite=True)

                    # Add to table
                    isot_time = masterflats_header['DATE']
                    time_dataframe = Time(isot_time, format='isot', scale='utc')
                    if element_list is None:
                        element_list = pd.DataFrame([[naxis1_expected, naxis2_expected,
                                                      x1, x2, y1, y2,
                                                      ccdbinx, ccdbiny,
                                                      int(name_dictionary[filtro]),
                                                      free_grisma, grisma_number,
                                                      flat_file_name, nights,
                                                      time_dataframe, time_dataframe.jd]],
                                                    columns=['Naxis1', 'Naxis2',
                                                             'x1', 'x2', 'y1', 'y2',
                                                             'Binning1', 'Binning2',
                                                             'filtro',
                                                             'free_grisma', 'num_grisma',
                                                             'flat_file_name', 'noche',
                                                             'tiempo_astropy', 'julian'])
                    else:
                        elemento_lista_ = pd.DataFrame([[naxis1_expected, naxis2_expected,
                                                         x1, x2, y1, y2,
                                                         ccdbinx, ccdbiny,
                                                         int(name_dictionary[filtro]),
                                                         free_grisma, grisma_number,
                                                         flat_file_name, nights,
                                                         time_dataframe, time_dataframe.jd]],
                                                       columns=['Naxis1', 'Naxis2',
                                                                'x1', 'x2', 'y1', 'y2',
                                                                'Binning1', 'Binning2',
                                                                'filtro',
                                                                'free_grisma', 'num_grisma',
                                                                'flat_file_name', 'noche',
                                                                'tiempo_astropy', 'julian'])
                        element_list = pd.concat([element_list, elemento_lista_], ignore_index=True)

            else:
                raise ValueError("Binning don't fit on both axis")

            if verbose >= 1:
                coord_lim = limites_imagen(*coordinates_image)
                imgdibujar(master_flats_colapsed, *coordinates_image, *coord_lim, verbose_=1)

            if interactive:
                input("Press Enter to continue...")

    return element_list


def make_master_flat(list_nights, list_bias, dir_lists, dir_data, dir_bias, dir_flats,
                     verbose, interactive, verbose_imagen=False):
    """
    Main function to make master flats.

    :param list_nights:
    :param list_bias:
    :param dir_lists:
    :param dir_data:
    :param dir_bias:
    :param dir_flats:
    :param verbose:
    :param interactive:
    :param verbose_imagen:
    :return:
    """

    i_night = 0
    df_flat = None
    for night in list_nights[:]:
        i_night += 1
        print('=== NIGHT ' + night + ' - (' + str(i_night) + '/' + str(len(list_nights)) + ') ===')

        lista_flats = read_list(dir_lists + night + '/' + 'LFlat.csv')

        sections, unique_section, section_counting, section_index, bin_sections, _ = create_unique_list(
            dir_data, night, lista_flats, header='CCDSEC', binning=True
        )

        sections_coordinates = np.zeros((len(unique_section), 4), dtype=int)
        for i in range(len(unique_section)):
            unique_coordinates = obtain_coordinates_ccd(unique_section[i])
            sections_coordinates[i, :] = [*unique_coordinates]

        df_flat_ = join_flat_images(night, unique_section, sections_coordinates, section_index,
                                    dir_bias, dir_data, dir_flats, lista_flats, list_nights, list_bias,
                                    verbose=verbose, interactive=interactive, verbose_images=verbose_imagen)

        if df_flat is None:
            df_flat = df_flat_
        else:
            df_flat = pd.concat([df_flat, df_flat_], ignore_index=True)

    return df_flat


def obtain_flats(dir_flats_, night_, list_nights, list_flats, x1, x2, y1, y2, b1, b2, id_filter_, max_search=10,
                 fill_flat=1, verbose=False):
    """
    Similarly to obtain_bias, searches the list of available flats for a flat with the same night the photo was taken.
    Also checks size nd binning. If it finds such a bias, it will use that one directly. If it does not find one, it
    will look for a different bias which fits the desired parameters for other nights. If there aren't any, it generates
    a flat bias.

    :param dir_flats_:
    :param night_:
    :param list_nights:
    :param list_flats:
    :param x1:
    :param x2:
    :param y1:
    :param y2:
    :param b1:
    :param b2:
    :param id_filter_:
    :param max_search:
    :param fill_flat:
    :param verbose:
    :return:
    """

    searched_flat_name = dir_flats_ + night_ + "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}-F{6:03d}.fits" \
        .format(x1, x2, y1, y2, b1, b2, id_filter_)
    exist = False
    take_another = False

    if searched_flat_name in list_flats:
        exist = True
    else:
        position = list_nights.index(night_)
        for i in range(1, max_search):
            if exist:
                break
            for mult in [-1, 1]:
                indx = i * mult
                position_new = position + indx
                if position_new >= len(list_nights):
                    # Reached end of the year. Doesn't search for next year. Can be fixed with a bigger folder and
                    # placing every night in a same folder
                    break
                night_ = list_nights[position_new]
                searched_flat_new = dir_flats_ + night_ + \
                    "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}-F{6:03d}.fits" \
                    .format(x1, x2, y1, y2, b1, b2, id_filter_)

                if searched_flat_new in list_flats:
                    searched_flat_name = searched_flat_new
                    exist = True
                    take_another = True
                    break

    if exist:
        if verbose:
            if take_another:
                print('Flat taken from another night. The taken flat is:')
                print(searched_flat_name)
            else:
                print('Flat exists.')
        searched_flat = fits.getdata(searched_flat_name, ext=0)
    else:
        if verbose:
            print('There are no nearby flats available. A filled bias has been generated.')

        # searched_flat = None
        # if searched_flat is None:
        #     raise ValueError('No hay definido un valor por defecto para el flat')

        naxis1_expected, naxis2_expected = obtain_naxis((x1, x2, y1, y2), b2, b1)

        searched_flat = np.full((naxis1_expected, naxis2_expected), fill_flat, dtype=float)

    return exist, searched_flat_name, searched_flat
