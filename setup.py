from setuptools import setup, find_packages

with open('requirements.txt') as requirements_file:
    install_requirements = requirements_file.read().splitlines()

# modify entry_points to use command line 
# {COMMAND NAME}={module path}:{function in the module}
setup(
    name="faersutil",
    version="0.0.1",
    description="a CLI package for handling FAERS data",
    author="tadahaya",
    packages=find_packages(),
    install_requires=install_requirements,
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "faersutil.preprocess=faersutil.preprocess:main",
            "faersutil.preprocess.parse=faersutil.preprocess:parse_xml",
            "faersutil.preprocess.clean=faersutil.preprocess:clean_and_merge",
            "faersutil.preprocess.curate=faersutil.preprocess:curate_drug",
            "faersutil.db=faersutil.database:main",
            "faersutil.db.drugdict=faersutil.database:update_drugdict",
            "faersutil.db.drug_rxn=faersutil.database:prep_drug_rxn",
            "faersutil.db.make_db=faersutil.database:make_database",
            "faersutil.calc=faersutil.calc:main",
        ]
    },
    classifiers=[
        'Programming Language :: Python :: 3.9',
    ]
)