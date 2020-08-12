"""Package installation"""

import setuptools

setuptools.setup(
    name='wrdists',
    entry_points={'console_scripts': ['wrdists = scripts.console_access:cli_main']},
    url='https://github.com/Gemma-Rate/wrdists',
    author='Gemma Rate',
    author_email='gemma.rate@gmail.com',
    install_requires=['numpy', 'astropy', 'matplotlib', 'pandas'],
    packages=setuptools.find_packages(),
    version='1.0',
    license='MIT',
    description='Calculates Galactic WR distances from Gaia DR2 parallaxes',
    long_description=open('README.md').read(),
)