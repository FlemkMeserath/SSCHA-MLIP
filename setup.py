from setuptools import setup, find_packages

VERSION = '0.3'
DESCRIPTION = 'MTP package for SSCHA'
LONG_DESCRIPTION = 'reworked a bit the cluster handling part and implemented the MTP package to include the SSCHA and MLIP'

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="utils_sscha_mlip",
        version=VERSION,
        author="Francesco Belli",
        author_email="<bellifra93@gmail.com>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=["python-sscha","CellConstructor"], # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'

        keywords=['python', 'SSCHA_MLIP'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Programming Language :: Python :: 3",
        ]
)

