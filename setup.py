from setuptools import setup
from setuptools import find_packages

def setup_package():

    from ultron import __version__
    REQUIRES = ['numpy', 'astropy>=1.0', 'pandas', 'matplotlib']
    META_DATA = dict(
        name='ultron',
        version=__version__,
        description='Image reduction for CAFOS instrument in CAHA',
        author='Enrique Galceran',
        author_email='egalcera@ucm.es',
        packages=find_packages('.'),
        entry_points={
            'console_scripts': [
                'ultron = ultron.ULTRON:main',
		'ultron-version = ultron.version:main'
            ],
            },
        setup_requires=['ultron'],
        install_requires=REQUIRES,
        zip_safe=False
        )

    setup(**META_DATA)


if __name__ == '__main__':
    setup_package()
