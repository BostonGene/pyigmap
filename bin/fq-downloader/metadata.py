from dataclasses import dataclass
from typing import List, Optional, Dict
from enum import Enum


class FileSource(str, Enum):
    FTP = "ftp"
    AWS = "aws"
    GCP = "gcp"
    NCBI = "ncbi"


@dataclass
class FileEntry:
    accession: str
    filename: str
    filetype: str
    filesize: Optional[int]
    filenumber: int
    md5: Optional[str]
    urltype: str
    url: str

    @classmethod
    def from_dict(cls, data: dict) -> "FileEntry":
        return cls(**data)


@dataclass
class Run:
    accession: str
    experiment: str
    study: str
    sample: str
    title: str
    attributes: Dict[str, str]
    files: Dict[FileSource, List[FileEntry]]

    @classmethod
    def from_dict(cls, data: dict) -> "Run":
        return cls(
            accession=data["accession"],
            experiment=data["experiment"],
            study=data["study"],
            sample=data["sample"],
            title=data["title"],
            attributes=data["attributes"],
            files=cls._parse_files(data["files"]),
        )

    @staticmethod
    def _parse_files(data: dict) -> Dict[FileSource, List[FileEntry]]:
        files: Dict[FileSource, List[FileEntry]] = {}
        for source in FileSource:
            entries = data.get(source, [])
            files[source] = [FileEntry.from_dict(f) for f in entries]
        return files


@dataclass
class Experiment:
    accession: str
    title: str
    platform: str
    instrument: str
    runs: Dict[str, Run]

    @classmethod
    def from_dict(cls, data: dict) -> "Experiment":
        runs = {k: Run.from_dict(v) for k, v in data["runs"].items()}
        return cls(
            accession=data["accession"],
            title=data["title"],
            platform=data["platform"],
            instrument=data["instrument"],
            runs=runs,
        )