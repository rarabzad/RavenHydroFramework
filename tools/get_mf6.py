#!/usr/bin/env python3
"""Download the MODFLOW 6 shared library (libmf6) that Raven's groundwater coupling loads at run time.

    python tools/get_mf6.py                  # the tested version (6.8.1) for this computer, into lib/mf6/
    python tools/get_mf6.py --latest         # the newest MODFLOW 6 release instead
    python tools/get_mf6.py --version 6.7.0  # a given release
    python tools/get_mf6.py --zip mf6.8.1_win64.zip   # use a release zip downloaded by hand (no internet needed)
    python tools/get_mf6.py --dest DIR       # put the library in DIR instead of lib/mf6/
    python tools/get_mf6.py --force          # replace a library that is already there

The release zips come from https://github.com/MODFLOW-ORG/modflow6/releases (for 6.8.1: mf6.8.1_linux.zip, mf6.8.1_win64.zip,
mf6.8.1_macarm.zip; mf6.X.Y_mac.zip for Intel Macs where released). Only the library is kept (libmf6.so, libmf6.dll or libmf6.dylib), with a
small text file naming its version and source. Raven finds it in lib/mf6/ of the source folder, or next to the Raven
executable, without any setting (see the manual, "Getting MODFLOW 6").

Standard Python 3 only; no other packages. MODFLOW 6 is public-domain software of the U.S. Geological Survey.
"""
import argparse
import hashlib
import json
import os
import platform
import re
import shutil
import sys
import tempfile
import urllib.request
import zipfile

TESTED = '6.8.1'    # the version the coupling was verified with (docs/VERIFICATION_REPORT_MF6.8.1.md)
REPO = 'MODFLOW-ORG/modflow6'
HERE = os.path.dirname(os.path.abspath(__file__))
DEST = os.path.join(os.path.dirname(HERE), 'lib', 'mf6')


def platform_info():
    """(release zip suffixes to try, library file name) for this computer"""
    s, m = platform.system(), platform.machine().lower()
    if s == 'Windows':
        if not (sys.maxsize > 2 ** 32 or m.endswith('64')):
            sys.exit('MODFLOW 6 is released for 64-bit Windows only')
        return ['win64'], 'libmf6.dll'
    if s == 'Darwin':
        if m in ('arm64', 'aarch64'):
            return ['macarm'], 'libmf6.dylib'
        return ['mac'], 'libmf6.dylib'   # Intel Macs: not built for every release
    if s == 'Linux':
        if m not in ('x86_64', 'amd64'):
            sys.exit('MODFLOW 6 releases a Linux library for x86_64 only (this is %s); build libmf6 from source: '
                     'https://github.com/MODFLOW-ORG/modflow6' % m)
        return ['linux'], 'libmf6.so'
    sys.exit('unsupported system %s: download libmf6 from https://github.com/%s/releases' % (s, REPO))


def request(url):
    return urllib.request.Request(url, headers={'User-Agent': 'Raven-get_mf6', 'Accept': 'application/octet-stream'})


def latest_version():
    """tag of the newest release: the GitHub API, else the redirect of /releases/latest"""
    try:
        with urllib.request.urlopen(request('https://api.github.com/repos/%s/releases/latest' % REPO), timeout=30) as r:
            return json.load(r)['tag_name'].lstrip('v')
    except Exception:
        pass
    try:
        with urllib.request.urlopen(request('https://github.com/%s/releases/latest' % REPO), timeout=30) as r:
            m = re.search(r'/releases/tag/v?([0-9][^/?#]*)', r.geturl())
    except Exception as e:
        sys.exit('could not reach GitHub to find the latest MODFLOW 6 release (%s); give --version' % e)
    if not m:
        sys.exit('could not find the latest MODFLOW 6 release; give --version')
    return m.group(1)


def download(url, path):
    """download url to path with a progress line; False if the file does not exist (404)"""
    try:
        with urllib.request.urlopen(request(url), timeout=60) as r, open(path, 'wb') as f:
            total = int(r.headers.get('Content-Length') or 0); got = 0
            while True:
                b = r.read(1 << 20)
                if not b:
                    break
                f.write(b); got += len(b)
                if total:
                    sys.stdout.write('\r  %s: %.0f of %.0f MB' % (os.path.basename(path), got / 1e6, total / 1e6)); sys.stdout.flush()
        if total:
            print()
        return True
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return False
        raise


def extract_library(zpath, libname):
    """(member name, bytes) of the library in a release zip"""
    with zipfile.ZipFile(zpath) as z:
        names = [n for n in z.namelist() if n.replace('\\', '/').split('/')[-1] == libname]
        if not names:
            sys.exit('%s holds no %s: is it a MODFLOW 6 release zip for this system?' % (zpath, libname))
        names.sort(key=lambda n: ('/bin/' not in n.replace('\\', '/'), len(n)))
        return names[0], z.read(names[0])


def main():
    ap = argparse.ArgumentParser(description='Download the MODFLOW 6 library for Raven (default: version %s).' % TESTED)
    g = ap.add_mutually_exclusive_group()
    g.add_argument('--version', help='MODFLOW 6 release, e.g. 6.8.1 (default: %s, the tested version)' % TESTED)
    g.add_argument('--latest', action='store_true', help='the newest MODFLOW 6 release')
    g.add_argument('--zip', help='a MODFLOW 6 release zip already downloaded (works offline)')
    ap.add_argument('--dest', default=DEST, help='folder for the library (default: lib/mf6 of the source folder)')
    ap.add_argument('--force', action='store_true', help='replace a library that is already there')
    a = ap.parse_args()

    suffixes, libname = platform_info()
    os.makedirs(a.dest, exist_ok=True)
    target = os.path.join(a.dest, libname)
    info = os.path.join(a.dest, 'libmf6_version.txt')
    if a.zip:
        version = (re.search(r'mf([0-9]+\.[0-9]+\.[0-9]+)', os.path.basename(a.zip)) or [None, 'unknown'])[1]
    else:
        version = latest_version() if a.latest else (a.version or TESTED)

    if os.path.exists(target) and not a.force:
        have = open(info).read().split('\n')[0].strip() if os.path.exists(info) else 'unknown version'
        if a.zip or version in have:
            print('%s is already there (%s); nothing to do (use --force to replace it)' % (target, have))
            return 0
        print('%s is already there (%s); use --force to replace it with MODFLOW %s' % (target, have, version))
        return 0

    tmp = tempfile.mkdtemp(prefix='get_mf6_')
    try:
        if a.zip:
            zpath, url = os.path.abspath(a.zip), a.zip
        else:
            zpath = None
            for suf in suffixes:
                url = 'https://github.com/%s/releases/download/%s/mf%s_%s.zip' % (REPO, version, version, suf)   # e.g. mf6.8.1_linux.zip
                p = os.path.join(tmp, os.path.basename(url))
                print('downloading %s' % url)
                if download(url, p):
                    zpath = p; break
            if zpath is None:
                extra = (' MODFLOW 6.8 and later are not built for Intel Macs; use an Apple Silicon Mac, an earlier '
                         'release (--version 6.6.3) or build libmf6 from source.' if suffixes == ['mac'] else '')
                sys.exit('MODFLOW %s has no release zip for this system (%s).%s See https://github.com/%s/releases'
                         % (version, ', '.join(suffixes), extra, REPO))
        member, data = extract_library(zpath, libname)
        part = target + '.part'
        with open(part, 'wb') as f:
            f.write(data)
        os.replace(part, target)                       # never leave a half-written library behind
        if platform.system() != 'Windows':
            os.chmod(target, 0o755)
        sha = hashlib.sha256(data).hexdigest()
        with open(info, 'w') as f:
            f.write('MODFLOW %s\nsource: %s (%s)\nsha256: %s\n' % (version, url, member, sha))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print('MODFLOW %s library: %s (%.1f MB)' % (version, target, len(data) / 1e6))
    if os.path.abspath(a.dest) != os.path.abspath(DEST):
        print('Raven finds it only next to its executable or in lib/mf6 of the source folder; otherwise set '
              'RAVEN_MF6_LIB=%s or :MF6Library in the .rvg file' % target)
    return 0


if __name__ == '__main__':
    sys.exit(main())
