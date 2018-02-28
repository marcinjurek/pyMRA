from setuptools import setup
setup(
    name = 'pyMRA',
    packages = ['pyMRA', 'pyMRA.data'],
    version = '0.7.7',
    description = 'Multi-resolution approximation for spatial Gaussian processes',
    author = 'Marcin Jurek',
    author_email = 'marcinjurek1988@gmail.com',
    url = 'https://github.com/marcinjurek/pyMRA',
    include_package_data=True,
    install_requires = ['numpy', 'matplotlib', 'scipy', 'numpy_indexed', 'sklearn'],
    #download_url = 'https://github.com/marcinjurek/pyMRA/archive/0.1.tar.gz',
    keywords = ['Gaussian process', 'multi-resolution', 'sparse', 'approximation', 'spatial', 'statistics'],
    classifiers = [],
)
