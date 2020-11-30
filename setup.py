import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pharmvip-guideline",
    version="0.0.1",
    author="Chanathip Sukritha",
    author_email="s.chanathip16@gmail.com",
    description="PharmVIP(Pharmacogenomic Variant Analysis and Interpretation Platform) - Guideline Module",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/cs16golf",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7.6',
    entry_points={
        "console_scripts":[
            "pharmvip_guideline=pharmvip_guideline.main:main"
        ]
    }
)
