"""Package installation"""

import setuptools

setuptools.setup(
    name='wrdists',
    url='https://github.com/Gemma-Rate/wrdists',
    author='Gemma Rate',
    author_email='gemma.rate@gmail.com',
    install_requires=['numpy', 'astropy', 'matplotlib', 'pandas'],
    packages=setuptools.find_packages(),
    entry_points ={'console_scripts': ['wrdists = wrdists.console_access:main']},
    version='1.1',
    license='MIT',
    description='Calculates Galactic WR distances from Gaia DR2 parallaxes',
    long_description=open('README.md').read(),
)