from setuptools import setup, find_packages

setup(
    name='qcmspycloud',
    version='0.2.0',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'streamlit',
        'streamlit-lottie>=0.0.5',
        'seaborn',
        "natsort"
        # other dependencies
    ],
    entry_points={
        'console_scripts': [
            'qcmspycloud=qcmspycloud.starter:starter',
        ],
    },
)
