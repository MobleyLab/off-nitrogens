"""
Example setup file
"""

import setuptools


if __name__ == "__main__":
    setuptools.setup(
        name='off_nitrogens',
        version="0.0.1",
        description='Determining improper dihedral parameters for trivalent nitrogens',
        author='Victoria Lim',
        author_email='limvt@uci.edu',
        url="https://github.com/MobleyLab/off_nitrogens",
        license='MIT',
        packages=setuptools.find_packages(),
        install_requires=[
            'numpy>=1.7',
        ],
        extras_require={
            'tests': [
                'pytest',
                'pytest-cov',
                'pytest-pep8',
                'tox',
            ],
        },

        tests_require=[
            'pytest',
            'pytest-cov',
            'pytest-pep8',
            'tox',
        ],

        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 3',
        ],
        zip_safe=True,
    )
