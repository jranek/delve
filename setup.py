from setuptools import setup

with open('README.md', 'r') as f:
    long_description = f.read()

setup(
    name='delve-fs',
    packages=['delve'],
    version='0.1.6',
    description='Feature selection for preserving biological trajectories from single-cell data',
    author='Jolene Ranek',
    author_email='ranekjs@gmail.com',
    url='https://github.com/jranek/delve',
    license='MIT',
    python_requires='3.9',
    long_description=long_description,
    long_description_content_type = 'text/markdown',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.9',
        'Operating System :: OS Independent',
    ],
    keywords=[
        'trajectory inference',
        'feature selection',
        'single-cell bioinformatics',
        'computational biology',
    ],
    install_requires=[
        'anndata>=0.7.6',
        'sketchKH==0.1.2',
        'numpy==1.26.4',
        'scipy>=1.7.1',
        'pandas==1.5.2',
        'scikit-learn>=0.23.2',
        'umap-learn==0.5.1',
        'matplotlib==3.7.1',
        'tqdm',
    ],
    ext_modules=[],
)