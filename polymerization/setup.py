from setuptools import setup, find_packages

setup(
    name='polymerization',             
    version='1.0.0',                 
    packages=find_packages(),      
    install_requires=[],           
    description='Gene fusions in populations',
    author='Jacob Krol',
    author_email='jacob.krol@colorado.edu',  
    classifiers=[                  
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
    ],
    entry_points={
    'console_scripts': [
        'poly-stix=polymerization.stix:main',
        # 'genefusion-stix_sharded=genefusion.stix_sharded:main',
        # 'genefusion-genefusion_stix_sharded=genefusion.genefusion_stix_sharded:main',
        # 'genefusion-genefusion_giggle=genefusion.genefusion_giggle:main',
        'poly-giggle_sharded=polymerization.giggle_sharded:main',
        'poly-stix_sharded=polymerization.stix_sharded:main',
        'poly-samplefusions=polymerization.samplefusions:main',
        'poly-index_els=polymerization.index_els:main',
        'poly-el2json=polymerization.el2json:main',
        'poly-g2ewdist=polymerization.g2ewdist:main',
        'poly-g2degst=polymerization.g2degst:main',
        'poly-fmt_giggle_fusion=polymerization.fmt_giggle_fusion:main'
    ],
    },
    python_requires='>=3.6',       
)
