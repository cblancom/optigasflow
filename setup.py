from setuptools import setup, find_packages

setup(name="Network",
    version="0.1", 
    description="Python script for compute gas natural flow", 
    author="cblancom",  
    author_email='cristian.blanco@utp.edu.co', 
    license="GPL",  
    url="https://github.com/cblancom/natural_gas_project",  
    packages=find_packages(),
    install_requires=[i.strip() for i in open("requirements.txt").readlines()],
)