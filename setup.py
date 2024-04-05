from setuptools import setup, find_packages

VERSION = '0.0.1'
DESCRIPTION = 'Compute gas natural flow'

# Setting up
setup(
    name="NGP",
    version=VERSION,
    author="cblancom",
    author_email="<cristian.blanco@utp.edu.co>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    package_dir={"": "optigasflow/src"},
    packages=find_packages(where="optigasflow/src"),
    install_requires=['cycler', 
                      'et-xmlfile',
                      'dill',
                      'fonttools', 
                      'gekko', 
                      'kiwisolver',
                      'matplotlib', 
                      'networkx',
                      'numpy', 
                      'openpyxl', 
                      'packaging', 
                      'pandas', 
                      'Pillow==9.3.0', 
                      'pyparsing',
                      'python-dateutil', 
                      'pytz', 
                      'scipy' , 
                      'six'],
    
    keywords=['python', 'natural gas'],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Programming Language :: Python :: 3",
    ]
)
