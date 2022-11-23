from setuptools import setup, find_packages
import codecs
import os

VERSION = '0.0.1'
DESCRIPTION = 'compute gas natural flow'

# Setting up
setup(
    name="gnf",
    version=VERSION,
    author="cblancom",
    author_email="<cristian.blanco@utp.edu.co>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    install_requires=['cycler', 'et-xmlfile', 'fonttools', 
                        'gekko', 'kiwisolver' ,'matplotlib', 
                        'networkx', 'numpy', 'openpyxl', 
                        'packaging', 'pandas', 'Pillow==9.3.0', 
                        'pyparsing','python-dateutil', 'pytz', 'scipy' , 'six'],
    keywords=['python', 'natural gas'],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Programming Language :: Python :: 3",
    ]
)