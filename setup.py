from setuptools import setup

setup(
    name='hmp2',
    version='0.0.1',
    description='Workflows, pipelines, and scripts for hmp2',
    packages=['hmp2'],
    install_requires=[
        'anadama', 'anadama_workflows', 'osdf-python', 'cutlass', 'toolz'
    ],
    zip_safe=False,
    classifiers=[
        "Development Status :: 2 - Pre-Alpha"
    ],
    entry_points={
        "anadama.pipeline": [
            ".internal_metadata = hmp2.internal_metadata:BroadMetadataPipeline",
            ".external_metadata = hmp2.external_metadata:ExternalMetadataPipeline",
        ]
    }
)
