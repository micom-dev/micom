"""Build a database of organism metabolic models."""

from urllib.parse import urlparse
import httpx
from rich.progress import (
    Progress,
    BarColumn,
    DownloadColumn,
    TextColumn,
    TransferSpeedColumn,
    TimeRemainingColumn,
)
import pandas as pd
from pathlib import Path
import os
from os import path
from zipfile import ZipFile
from .constants import DB_URL, MEDIA_URL


def load_manifest(folder):
    """Get the manifest from a model DB."""
    mpath = path.join(folder, "manifest.csv")
    if not path.exists(mpath):
        raise ValueError(
            "No manifest found. `%s` does not look like a valid "
            "model database." % folder
        )
    manifest = pd.read_csv(path.join(folder, "manifest.csv"))
    if "file" not in manifest.columns:
        raise ValueError("Invalid manifest for model database :(")
    manifest.file = [path.join(folder, f) for f in manifest.file]
    return manifest


def load_zip_model_db(artifact, extract_path):
    """Prepare a model database for use."""
    if not path.exists(extract_path):
        os.mkdir(extract_path)
    with ZipFile(artifact) as zf:
        zf.extractall(extract_path)
    manifest = load_manifest(extract_path)
    manifest["file"] = [path.join(extract_path, f) for f in manifest.file]
    return manifest


def get_database(url: str, out: Path, what: str = "taxa") -> Path:
    """Get a database from several locations.

    If the database is a local file it will be used directly. If it is a URL it will be downloaded to the specified location.

    Parameters
    ----------
    url : str
        The URL of the database to download or a locally downloaded database.
    out : Path
        The path to the folder where the database should be downloaded.
    what : str
        The type of database to download.

    Returns
    -------
    Path
        The path to the downloaded database.

    """

    if url is None:
        return None

    base = DB_URL if what == "taxa" else MEDIA_URL

    progress = Progress(
        TextColumn("{task.fields[database]}", justify="right"),
        BarColumn(bar_width=None),
        "[progress.percentage]{task.percentage:>3.1f}%",
        DownloadColumn(),
        TransferSpeedColumn(),
        TimeRemainingColumn(),
    )

    up = urlparse(url)
    if up.scheme == "default" and up.netloc:
        dl = base % up.netloc
        loc = out / up.netloc
    elif (up.scheme in ["http", "https"]) and up.netloc:
        dl = url
        loc = out / Path(up.path).name
    elif up.scheme == "file":
        return Path(up.netloc)
    else:
        return Path(url)

    out.mkdir(parents=True, exist_ok=True)

    with (
        progress,
        httpx.stream(
            method="GET",
            url=dl,
            follow_redirects=True,
            timeout=60,
        ) as response,
        open(loc, "wb") as data,
    ):
        response.raise_for_status()
        print(f"Connected. Downloading model database to {loc}.")
        task_id = progress.add_task(
            description="download model database",
            database=up.netloc,
            total=int(response.headers.get("Content-Length", None)),
        )
        for chunk in response.iter_bytes():
            data.write(chunk)
            progress.update(task_id=task_id, advance=len(chunk))

        return loc
