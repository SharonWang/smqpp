from setuptools import setup, find_packages

setup(
	name='smqpp',
	version='0.1.5',
	description='Smartseq2 preprocessing toolkit',
	url='https://github.com/SharonWang/smqpp',
	packages=find_packages(exclude=['docs', 'figures', 'examples']),
	install_requires=['numpy', 'matplotlib', 'pandas', 'anndata', 'scipy', 'statsmodels', 'patsy', 'umap'],
	author='Xiaonan Wang',
	author_email='xw251@cam.ac.uk',
	license='MIT'
)
