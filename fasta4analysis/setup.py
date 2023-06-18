from setuptools import setup, find_packages

with open('requirements.txt') as f:
	requirements = f.readlines()

long_description = 'DNA Sequences Data Analysis tool using fasta, txt, and csv \
files data to create histogram, barh and pie chart plots'

setup(
		name ='fasta4analysis',
		version ='1.0.0',
		author ='Swati Mishra',
		author_email ='mishraswati1387@gmail.com',
		url ='https://github.com/GitCodeSM',
		description ='Package for fasta matplotlib Data Analysis',
		long_description = long_description,
		long_description_content_type ="text/markdown",
		license ='MIT', 
		packages = find_packages(),
		entry_points ={
			'console_scripts': [
				'fastacli = fasta4analysis.fastacli:main'
			]
		},
		classifiers =(
			"Programming Language :: Python :: 3.10.4",
			"License :: OSI Approved :: MIT License",
			"Operating System :: OS Independent",
		),
		keywords ='fasta matplotlib python package',
		install_requires = requirements,
		zip_safe = False
)
