"""
step00.py: for base implementation
"""
import hashlib
import hmac
import os
import pickle
import re
import tarfile
import tempfile
import typing

secret = bytes("asdf", "UTF-8")
tmpfs = "/tmpfs"


def file_list(path: str) -> typing.List[str]:
    """
    file_list: return a list of files in path
    """
    return list(filter(lambda x: os.path.isfile(x), list(map(lambda x: os.path.join(path, x), os.listdir(path)))))


def directory_list(path: str) -> typing.List[str]:
    """
    directory_list: return a list of directories in path
    """
    return list(filter(lambda x: os.path.isdir(x), list(map(lambda x: os.path.join(path, x), os.listdir(path)))))


def make_pickle(path: str, data: typing.Any) -> None:
    """
    make_pickle: create a pickle
    """
    if not path.endswith(".tar.xz"):
        raise ValueError("Path must end with .tar.xz")

    pkl = pickle.dumps(data, protocol=pickle.HIGHEST_PROTOCOL)
    key = hmac.new(secret, pkl, hashlib.sha512).digest()

    with tempfile.TemporaryDirectory(dir=tmpfs) as tmp_dir:
        with open(os.path.join(tmp_dir, "data.pkl"), "wb") as f:
            f.write(pkl)
        with open(os.path.join(tmp_dir, "key.txt"), "wb") as f:
            f.write(key)

        with tarfile.open(path, "w:xz") as tar:
            tar.add(os.path.join(tmp_dir, "data.pkl"), arcname="data.pkl")
            tar.add(os.path.join(tmp_dir, "key.txt"), arcname="key.txt")


def read_pickle(path: str) -> typing.Any:
    """
    read_pickle: read a pickle file
    """
    if not path.endswith(".tar.xz"):
        raise ValueError("Path must end with .tar.xz")
    if not tarfile.is_tarfile(path):
        raise ValueError("Path cannot be read as a tar file")

    with tempfile.TemporaryDirectory(dir=tmpfs) as tmp_dir:
        with tarfile.open(path, "r:xz") as tar:
            tar.extractall(tmp_dir)

        with open(os.path.join(tmp_dir, "data.pkl"), "rb") as f:
            pkl = f.read()
        with open(os.path.join(tmp_dir, "key.txt"), "rb") as f:
            key = f.read()

    if not hmac.compare_digest(hmac.new(secret, pkl, hashlib.sha512).digest(), key):
        raise ValueError("Data is not valid")

    return pickle.loads(pkl)


def get_patient(ID: str) -> str:
    return re.findall(r"(^(cn)?\d+)", ID)[0][0]
