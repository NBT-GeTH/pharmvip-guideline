import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pharmvip-guideline",
    version="0.0.5",
    author="Chadapohn Chaosrikul, Krittin Phornsiricharoenphant, Chanathip Sukritha",
	author_email="chadapohn.chaosrikul@gmail.com, oatkrittin@gmail.com, s.chanathip16@gmail.com",
	description="PharmVIP(Pharmacogenomic Variant Analysis and Interpretation Platform) - Guideline Module",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7.6',
    entry_points={
        "console_scripts":[
            "pharmvip_guideline=pharmvip_guideline.main:main"
        ]
    }
)
