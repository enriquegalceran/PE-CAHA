import csv


def tuple2coordinates(tupla):
    return tupla[0], tupla[1], tupla[2], tupla[3]


def save_file_csv(csvfile, res):
    with open(csvfile, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        for val in res:
            writer.writerow([val])
