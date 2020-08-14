"""Package installation"""

import setuptools

setuptools.setup(
    name='wrdists',
    scripts=['wrdists/scripts/console_access'],
    url='https://github.com/Gemma-Rate/wrdists',
    author='Gemma Rate',
    author_email='gemma.rate@gmail.com',
    install_requires=['numpy', 'astropy', 'matplotlib', 'pandas'],
    packages=['wrdists','wrdists.test'],
    version='1.0',
    license='MIT',
    description='Calculates Galactic WR distances from Gaia DR2 parallaxes',
    long_description=open('README.md').read(),
)