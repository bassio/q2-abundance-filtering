from setuptools import setup, find_packages


setup(
    name="q2-abundance-filtering",
    version="0.01",
    packages=find_packages(),
    author="Ahmed Bassiouni",
    author_email="ahmedbassi@gmail.com",
    description="QIIME 2 plugin for abundance-filtering of sequences according to the method of Wang et al.",
    license='BSD-3-Clause',
    url="https://github.com/bassio/q2-abundance-filtering",
    entry_points={
        "qiime2.plugins":
        ["q2-abundance-filtering=q2_abundance_filtering.plugin_setup:plugin"]
    },
    package_data={
        'q2_abundance_filtering': ['citations.bib']
        },
    zip_safe=False,
) 
