from setuptools import setup, find_packages

setup(
    name='genefusion',             
    version='0.1',                 
    packages=find_packages(),      
    install_requires=[],           
    description='Gene fusion helpers',
    author='Jacob Krol',
    author_email='jacob.krol@colorado.edu',  
    classifiers=[                  
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
    ],
    entry_points={
    'console_scripts': [
        'genefusion-stix=genefusion.stix:main',
        'genefusion-stix_sharded=genefusion.stix_sharded:main',
        'genefusion-genefusion_stix_sharded=genefusion.genefusion_stix_sharded:main',
        'genefusion-genefusion_giggle=genefusion.genefusion_giggle:main',
    ],
    },
    python_requires='>=3.6',       
)
