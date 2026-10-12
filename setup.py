import re
from pathlib import Path

from setuptools import setup
from setuptools_rust import Binding, RustExtension

version = re.search(r'^version = "(.+)"', Path("Cargo.toml").read_text(), re.MULTILINE)[
    1
]

setup(
    version=version,
    packages=["groebner"],
    package_dir={"": "python"},
    rust_extensions=[
        RustExtension(
            "groebner.groebner",
            binding=Binding.NoBinding,
            features=["ffi"],
            py_limited_api=True,
        )
    ],
    options={"bdist_wheel": {"py_limited_api": "cp39"}},
    zip_safe=False,
)
