"""Shared local release-file enumeration. No scientific calculations or network."""
from __future__ import annotations
import hashlib
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
EXCLUDED_DIRS = {'.git', '.venv', 'venv', '__pycache__', '.pytest_cache', '.matplotlib-cache'}
EXCLUDED_SUFFIXES = {'.pyc', '.pyo', '.aux', '.out', '.toc', '.tmp'}

def payload_files(root: Path = ROOT) -> list[Path]:
    files = []
    for path in sorted(root.rglob('*')):
        relative = path.relative_to(root)
        if not path.is_file() or path.is_symlink():
            continue
        if relative.parts[0] in {'results', 'dist'} or any(p in EXCLUDED_DIRS for p in relative.parts):
            continue
        if path.suffix.lower() in EXCLUDED_SUFFIXES or path.name in {'.DS_Store', 'Thumbs.db'}:
            continue
        files.append(path)
    return files

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()

def write_checksums(root: Path = ROOT) -> Path:
    lines = [sha256(path) + '  ' + path.relative_to(root).as_posix()
             for path in payload_files(root) if path.name != 'SHA256SUMS']
    result = root/'SHA256SUMS'
    result.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return result

def check_checksums(root: Path = ROOT) -> list[str]:
    manifest = root/'SHA256SUMS'
    if not manifest.is_file():
        return ['SHA256SUMS missing']
    errors, listed = [], set()
    for line in manifest.read_text().splitlines():
        expected, name = line.split('  ', 1)
        listed.add(name)
        path = (root/name).resolve()
        if root.resolve() not in path.parents or not path.is_file():
            errors.append('Missing or unsafe file: ' + name)
        elif sha256(path) != expected:
            errors.append('Checksum mismatch: ' + name)
    actual = {p.relative_to(root).as_posix() for p in payload_files(root) if p.name != 'SHA256SUMS'}
    errors.extend('Unlisted file: ' + name for name in sorted(actual-listed))
    errors.extend('Listed but excluded/missing file: ' + name for name in sorted(listed-actual))
    return errors
