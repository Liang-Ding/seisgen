import setuptools

setuptools.setup(
    name="seisgen",
    version="0.1.7",
    author="Liang Ding",
    author_email="myliang.ding@mail.utoronto.ca",
    description="The theoretical seismic waveform generation tool",
    long_description="The theoretical seismic waveform generation tool based on the stored Strain Green's Tensor (SGT)"
                     " database that is built by using the 3D background model and SPECFEM3D package.",
    long_description_content_type="text/markdown",
    url="https://github.com/Liang-Ding/seisgen",
    project_urls={
        "Bug Tracker": "https://github.com/Liang-Ding/seisgen/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    keywords=[
        "seismology"
    ],
    # package_dir={"": "seisgen"},
    python_requires='>=3.6.0',
    install_requires=[
        "numpy", "scipy", "obspy",
        "h5py",
    ],
    # packages=setuptools.find_packages(where="seisgen"),
    packages=setuptools.find_packages(),
)
