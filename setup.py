import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PBLMM",
    version="2.1.2",
    author="Kevin Klann",
    author_email="klann@em.uni-frankfurt.de",
    # Updated by = "Süleyman Bozkurt",
    # Updated date = '25/11/2025'
    description="Python package for linear mixed model for peptide level of proteomics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/science64/PBLMM",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)