"""Getting and finding the MODFLOW 6 library without any setting.

1. tools/get_mf6.py --zip: a release zip (built here from the configured library, laid out like the MODFLOW release
   zips: mf6.X.Y_<platform>/bin/libmf6.*) is unpacked into a folder, with its version file; a second call does nothing.
2. Raven with RAVEN_MF6_LIB unset finds the library next to its executable, and in lib/mf6 one folder up (the layout of
   src/Raven.exe or build/Raven in a source checkout); GWModelSummary.txt names the library it used.
3. With no library anywhere, the run stops with a message that lists the places searched and names tools/get_mf6.py.
"""
import os, sys, shutil, subprocess, zipfile, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import testconfig as tc

tc.need_mf6()
W = tc.work('mf6lib')
shutil.rmtree(W, ignore_errors=True); os.makedirs(W)
name = os.path.basename(tc.MF6LIB) if os.path.basename(tc.MF6LIB).startswith('libmf6.') else \
    {'win32': 'libmf6.dll', 'darwin': 'libmf6.dylib'}.get(sys.platform, 'libmf6.so')
ok = True


def check(label, cond, detail=''):
    global ok
    print('%-62s %s %s' % (label, 'ok' if cond else 'FAIL', detail)); ok &= bool(cond)


# 1. offline install from a release zip
plat = {'win32': 'win64', 'darwin': 'macarm'}.get(sys.platform, 'linux')
z = os.path.join(W, 'mf6.8.1_%s.zip' % plat)
with zipfile.ZipFile(z, 'w', zipfile.ZIP_DEFLATED) as f:
    f.write(tc.MF6LIB, 'mf6.8.1_%s/bin/%s' % (plat, name))
    f.writestr('mf6.8.1_%s/doc/readme.txt' % plat, 'not the library')
dest = os.path.join(W, 'lib', 'mf6')
get = [sys.executable, os.path.join(tc.REPO, 'tools', 'get_mf6.py'), '--zip', z, '--dest', dest]
r = subprocess.run(get, capture_output=True, text=True)
lib = os.path.join(dest, name)
same = os.path.exists(lib) and open(lib, 'rb').read() == open(tc.MF6LIB, 'rb').read()
check('get_mf6.py --zip unpacks the library', r.returncode == 0 and same, r.stderr.strip()[-200:])
ver = open(os.path.join(dest, 'libmf6_version.txt')).read() if os.path.exists(os.path.join(dest, 'libmf6_version.txt')) else ''
check('  version file names MODFLOW 6.8.1 and a sha256', 'MODFLOW 6.8.1' in ver and 'sha256:' in ver)
r = subprocess.run(get, capture_output=True, text=True)
check('  a second call leaves it alone', r.returncode == 0 and 'nothing to do' in r.stdout)


# 2./3. Raven finds the library by itself
def case(tag, exe_dir):
    d = os.path.join(W, tag)
    tc.copy_case(tc.LIARD, d)
    t = open(d + '/Liard.rvi').read(); t = re.sub(r':Duration\s+\S+', ':Duration 2', t); open(d + '/Liard.rvi', 'w').write(t)
    os.makedirs(exe_dir, exist_ok=True)
    exe = os.path.join(exe_dir, os.path.basename(tc.RAVEN)); shutil.copy2(tc.RAVEN, exe)
    env = {k: v for k, v in tc.ENV.items() if k != 'RAVEN_MF6_LIB'}
    os.makedirs(d + '/out', exist_ok=True)
    r = subprocess.run([exe, 'Liard', '-o', 'out/'], cwd=d, capture_output=True, text=True, env=env, timeout=1800)
    msg = (re.findall(r'Exiting Gracefully: (.*)', r.stdout + r.stderr) or ['(no exit message)'])[-1]
    summ = open(d + '/out/GWModelSummary.txt').read() if os.path.exists(d + '/out/GWModelSummary.txt') else ''
    used = (re.findall(r'MODFLOW 6 library\s*: (.*)', summ) or [''])[0]
    return msg, used


tree = os.path.join(W, 'checkout')                       # checkout/src/Raven.exe + checkout/lib/mf6/libmf6.*
os.makedirs(os.path.join(tree, 'lib', 'mf6')); shutil.copy2(lib, os.path.join(tree, 'lib', 'mf6', name))
msg, used = case('in_lib_mf6', os.path.join(tree, 'src'))
check('found in lib/mf6 one folder up (src/Raven.exe)', 'Successful' in msg and 'lib/mf6' in used.replace('\\', '/'), used)

beside = os.path.join(W, 'beside')                       # the library next to the executable
os.makedirs(beside); shutil.copy2(lib, os.path.join(beside, name))
msg, used = case('beside_exe', beside)
check('found next to the executable', 'Successful' in msg and os.path.dirname(used.split(' (')[0]) == beside, used)

msg, used = case('missing', os.path.join(W, 'iso', 'x', 'y', 'alone'))   # no lib/mf6 within two folders up
check('no library: the message lists the places and get_mf6.py',
      ('get_mf6.py' in msg) and ('RAVEN_MF6_LIB' in msg) and ('lib/mf6' in msg.replace('\\', '/')), msg[:120])

print('library checks: %s' % ('ALL PASSED' if ok else 'FAILED'))
sys.exit(0 if ok else 1)
