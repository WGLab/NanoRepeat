from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()
    
setup(
    name='NanoRepeat',
    version='1.8.3',
    description='NanoRepeat: quantification of Short Tandem Repeats (STRs) from long-read sequencing data',
    long_description = readme,
    long_description_content_type = 'text/markdown',
    url='https://github.com/WGLab/NanoRepeat',
    author='Li Fang, Kai Wang',
    author_email='fangli9@mail.sysu.edu.cn',
    license='MIT',
    packages = find_packages("src"),
    package_dir = {"": "src"},
    data_files = [("", ["LICENSE"])],
    scripts=['src/NanoRepeat/nanoRepeat.py', 'src/NanoRepeat/nanoRepeat_joint.py'],
    install_requires=['matplotlib>=3.4.0',
                      'numpy>=1.21.6', 
                      'scikit-learn>=0.22.1',
                      'pysam>=0.17.0',
                      'pyminimap2>=2.30.0',
                      'levenshtein>=0.25']
)
