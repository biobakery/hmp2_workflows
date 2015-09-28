from setuptools import setup

setup(
    name='hmp2',
    version='0.0.1',
    description='Workflows, pipelines, and scripts for hmp2',
    packages=['hmp2'],
    zip_safe=False,
    install_requires=[
        'anadama',
        'anadama_workflows'
    ],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha"
    ]
)
