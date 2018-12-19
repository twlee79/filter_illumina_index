import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="filter_illumina_index",
    version="1.0.1",
    author="Tet Woo Lee",
    author_email="developer@twlee.nz",
    description="Filter a Illumina FASTQ file based on index sequence",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/twlee79/filter_illumina_index",
    packages=setuptools.find_packages(),
    install_requires=[
        'biopython',
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    entry_points={
        'console_scripts': [
            'filter_illumina_index = filter_illumina_index.filter_illumina_index:main',
        ],
    },
    data_files=[("", ["LICENSE"])],
)
