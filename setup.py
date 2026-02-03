import setuptools

setuptools.setup(
    name='cafea',
    version='0.0.0',
    description='Coffea analysis framework en Asturias',
    packages=setuptools.find_packages(),
    # Include data files (Note: "include_package_data=True" does not seem to work)
    package_data={
        "cafea" : [
            "cfg/*.cfg",
            "json/*",
            "data/scaleFactors/*.root",
            "data/fromTTH/fakerate/*.root",
            "data/fromTTH/fliprates/*.root",
            "data/fromTTH/lepSF/*/*/*.root",
            "data/fromTTH/lepSF/*/*/*/*.root",
            "data/JEC/*.txt",
            "data/btagSF/UL/*.pkl.gz",
            "data/btagSF/UL/*.csv",
            "data/btagSF/*.csv",
            "data/pileup/*.root",
            "data/goldenJsons/*.txt",
        ],
    }
)

