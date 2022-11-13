# import pandas as pd
# from network import *

# path_file = "/home/cristianblanco/projects/natural_gas_project/database/ng_case8.xlsx"
# data = pd.ExcelFile(path_file)
# p = Network(data) 

# print(p.get_values())


from setuptools import setup, find_packages


setup(
    name="foo",
    version="1.0",
    packages=find_packages(),
)