from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='Catactor',
    version='0.1.2',
    description="scATAC-seq analysis using meta-analytic marker genes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/carushi/Catactor/",
    author="carushi",
    author_email="rkawaguc@cshl.edu",
    license="MIT",
    packages=find_packages(include=['Catactor', 'Catactor.*']),
    package_data={'Catactor': ['./marker_genes/others/marker_name_list.csv', './marker_genes/others/rank_gene_list_celltype.csv']},
    include_package_data=True
)    

