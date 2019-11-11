# PE-CAHA

This proyect uses Python 3 to automate a pipeline and generate a series of calibrated telescope images for the CAFOS instrument on the 2.5m telescope in CAHA, Andalusia, Spain.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

You will need a python environment for version >=3 (preferrably 3.7) with several packages:

```
apt-get update
apt-get install python-virtualenv
virtualenv -p /usr/bin/python3 path/to/virtual/environment/folder/CAHA
```

To activate the environment:
```
source path/to/virtual/environment/folder/CAHA/bin/activate
```

Once activated, you should see something along these lines.
```
(CAHA) root@host2:~#
```

Once activated, you should install several python modules (Lista no comprobada que estén todos los módulos)

```
pip install numpy pandas scipy math astropy json csv
```

### Installing

To install this code, you will need to download the github repository and install it using the setup.py code to make it available in your machine.

Clone the github repository

```
$ git clone https://github.com/enriquegalceran/PE-CAHA.git
```

Go to the folder where you cloned the repository and install using python

```
cd path/to/repository
python setup.py make
python setup.py install
```
(You probably will need to add the directory where you installed the code to PATH)

## Running the tests

(Debería añadir una noche en el repositorio para que se pueda usar como prueba. Así, los valores por defecto para el argsparse pueden hacer los datos demo)

### Generating Bias files

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

# How to use ULTRON

ULTRON ()_Unidad de Limpieza y Tratamiento de Resultados Observacionales Nativo_) is an automated system, that given a folder with several subfolders for each night of observation will analyse each folder and generate a directory with the calibrated and reduced images. At the moment, ULTRON is made specifically for the CAFOS instrument in the 2.5m telescope in CAHA.

To operate ULTRON, you will need to add multiple directories into the command line while executing:

```
ultron -db path/to/masterbias/output -df path/to/masterflats/output -dl path/to/aux/list/directory \
              -de path/to/calibrated/output
```

(He pensado, para que el código estuviera completo, que en realidad con que le des una dirección para esto, sería suficiente: Le das la dirección de archivos de salida y si no están allí dentro las carpetas las genera; pero que de todas formas, si defines la carpeta sobreescriba esa salida automática)

While the code is being executed, auxiliary files (e.g., dataframes and lists that identify each image) will be stored. These can afterwards be deleted. These dataframes help in the later stages of the code by generating an easy-to-acess dataframe with the available calibration images (it's easier to sort a dataframe than opening and closing multiple .fits files, specially if that was has already been done).

Additionally, you can command the orquestrator to skip different modules (i.e., Masterbias or  Masterflat generation, or reducing the final images), as well as forcing to execute a speciffic module anew. This is usefull if you are just creating calibration images or if a calibration image has not been done correctly.

To do this, add the different arguments to skip or force:

```
ultron --nobias    # To skip masterbias
ultron --sibias    # To force masterbias
ultron --noflat    # To skip masterflas
ultron --siflat    # To force masterflat
ultron --noreducc  # To skip final reduction of the images
```

## Additional useful information

The clasification between different type of images is done automatically. As such, it can miss-clasify images. However, the classification the observer should use the tools available to try to minimize possible wrong labels. If there is a conflict between multiple types of images (for example, a bias labeled originally as a flat), the clasification module will try to label it correctly. In those cases, an error log will be generated for each night with a list of images where there was an ambiguity. I recommend that although in theory you can simply run the code and it should fix itself, to run first of all just the clasifier (using all the skip options) and check mannually the discrepancies (there should just be handfull of them. During the elaboration of this code, from 6440 test images, under 30 where wrongly classified, and most of them where from a single night where the observer forgot to switch from bias to science)).

If you want to generate anew the lists for a specific night, the best method is to delete the FOLDER with all the lists, not just the content in it. The next time the code starts, it will create the lists from scratch. Use this if you have fixed mannualy the headers for a faulty fits.

_Experimental:_ In theory, with the ```--calysci``` option, it should work with images observed with CAFOS which don't have the naming convention used in the Spanish Virtual Observatory. I have tried to use this option to calibrate the images which have been directly downloaded from CAHA and they worked, but there are still tests to be made.

_For future contributors:_ An initial Progressbar has been added (called ```progressbar.py```), which generates a progressbar in the last line in the terminal. The idea was to generate said line while the code was working and thus specify how long it is going to take. The idea is simple: instead of ending in a new line, with '\r' the cursor starts from the beginning, so it overwrites the last iteration of the progressbar. However, if there are other lines being printed, will will ocassionally bug and create multiple instances of the progressbar. It did work if there was an update to the progressbar after EVERY print: The progressbar is printed and when there is a new print it overwrites the arrow and if there is a new print of the progressbar, it will write a new line (by default) and the cursor will be again at the start of the line ready to be overwritten. This works if the computer is fast enough and the monitor has a fast enough refreshrate (on 50Hz monitors it twinkles and the effect is lost).

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
