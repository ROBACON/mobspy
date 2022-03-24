import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt") as fh:
    install_requires = fh.read().strip().split("\n")

with open("mobspy/_version.py") as fh:
    version = fh.readlines()[-1].split()[-1].strip("\"'")

setuptools.setup(
    name="mobspy",
    packages=["mobspy"],
    version=version,
    description="A Query-Based Language for Chemical Reaction Networks",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ROBACON/mobspy",
    python_requires=">=3.6",
    install_requires=install_requires
)
