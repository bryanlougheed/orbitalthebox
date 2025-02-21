from setuptools import setup, find_packages

setup(
    name='orbitalthebox',
    version='0.1.0',
    description='Calculating irradiance for Earth for various astronomical parameters.',
    author='Bryan C. Lougheed',
    author_email='bryan.lougheed@outlook.com',
    packages=find_packages(),
    install_requires=[
        # List your package dependencies here
        'numpy',
        'pandas',
        'os',
    ],
)
