# Keep the vendored IRI Fortran/data files in sync with irimodel.org before
# any iriXX extension_module() below references them -- see iri/sync_iri.py's
# own docstring for what this does and does not touch, and why. Failures
# (no network, site unreachable) are swallowed inside the script itself, so
# this never blocks a build -- it just leaves the last-synced (committed)
# copy in place; check: false here is only a second line of defense.
sync_iri_result = run_command(py3, files('iri/sync_iri.py'), check: false)
if sync_iri_result.returncode() != 0
  warning('iri/sync_iri.py failed: ' + sync_iri_result.stderr().strip())
endif

rocsat_module = custom_target('rocsat_module',
    output: 'rocsatfmodule.c',
    input: ['rocsat/rocsatf.pyf'],
    command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('rocsatf',
    ['rocsat/vfjmodroc1.F90',
     rocsat_module,
    ],
    include_directories: [inc_np, inc_f2py],
    dependencies: fortranobject_dep,
    install: true,
    link_language: 'fortran',
    subdir: 'mpyricalspace'
)

scherliessfejer_module = custom_target('scherliessfejer_module',
   output: 'scherliessfejerfmodule.c',
   input: ['scherliessfejer/scherliessfejerf.pyf'],
   command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('scherliessfejerf',
   ['scherliessfejer/StormDrift.F90',
    scherliessfejer_module,
   ],
   include_directories: [inc_np, inc_f2py],
   dependencies: fortranobject_dep,
   install: true,
   link_language: 'fortran',
   subdir: 'mpyricalspace'
)

# Keep hwm14/hwm07's own files in sync with map.nrl.navy.mil before either
# extension_module() below references them -- see hwm/sync_hwm.py's own
# docstring for what this does and does not touch, and why.
sync_hwm_result = run_command(py3, files('hwm/sync_hwm.py'), check: false)
if sync_hwm_result.returncode() != 0
  warning('hwm/sync_hwm.py failed: ' + sync_hwm_result.stderr().strip())
endif

hwm14_module = custom_target('hwm14_module',
   output: 'hwm14fmodule.c',
   input: ['hwm/hwm14/hwm14f.pyf'],
   command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('hwm14f',
   ['hwm/hwm14/hwm14.f90',
    'hwm/hwm14/hwm14_batch.f90',
    hwm14_module,
   ],
   include_directories: [inc_np, inc_f2py],
   dependencies: fortranobject_dep,
   install: true,
   link_language: 'fortran',
   subdir: 'mpyricalspace'
)

hwm07_module = custom_target('hwm07_module',
   output: 'hwm07fmodule.c',
   input: ['hwm/hwm07/hwm07f.pyf'],
   command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('hwm07f',
   ['hwm/hwm07/hwm07.01d.f90',
    'hwm/hwm07/hwm07_batch.f90',
    hwm07_module,
   ],
   include_directories: [inc_np, inc_f2py],
   dependencies: fortranobject_dep,
   install: true,
   link_language: 'fortran',
   subdir: 'mpyricalspace'
)

hwm93_module = custom_target('hwm93_module',
   output: 'hwm93fmodule.c',
   input: ['hwm/hwm93/hwm93f.pyf'],
   command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('hwm93f',
   ['hwm/hwm93/hwm93.f',
    'hwm/hwm93/hwm93_batch.f90',
    hwm93_module,
   ],
   include_directories: [inc_np, inc_f2py],
   dependencies: fortranobject_dep,
   install: true,
   link_language: 'fortran',
   subdir: 'mpyricalspace'
)

hltwim_module = custom_target('hltwim_module',
   output: 'hltwimfmodule.c',
   input: ['hltwim/hltwimf.pyf'],
   command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('hltwimf',
   ['hltwim/hltwim.f90',
    hltwim_module,
   ],
   include_directories: [inc_np, inc_f2py],
   dependencies: fortranobject_dep,
   install: true,
   link_language: 'fortran',
   subdir: 'mpyricalspace'
)

m_dep = cc.find_library('m', required : true)
gsl_dep = dependency('gsl')            # external: brew install gsl / apt install libgsl-dev / conda install gsl

py3.extension_module('jvdm1',
   ['jvdm1/jvdm1.c',
   'jvdm1/drift_basis.c',
   'jvdm1/drift_calc.c',
   'jvdm1/ndlinear.c',                 # P. Alken's gsl_multifit_ndlinear, vendored (was cc.find_library('ndlinear'))
   ],
   include_directories: include_directories('jvdm1'),
   dependencies: [m_dep, gsl_dep],
   install: true,
   subdir: 'mpyricalspace'
)

# Keep eejm1/eejm2's own files in sync with geomag.colorado.edu before either
# extension_module() below references them -- see eejm/sync_eejm.py's own
# docstring for what this does and does not touch, and why.
sync_eejm_result = run_command(py3, files('eejm/sync_eejm.py'), check: false)
if sync_eejm_result.returncode() != 0
  warning('eejm/sync_eejm.py failed: ' + sync_eejm_result.stderr().strip())
endif

# P. Alken's Equatorial Electrojet (v1, v2) and Equatorial Electric Field models.
# Each *.c wrapper (L. Navarro) batch-loops the calc routines behind one numpy call.
py3.extension_module('eejm1',
   ['eejm/eejm1/eejm1.c', 'eejm/eejm1/eej_basis.c', 'eejm/eejm1/eej_calc.c',
    'eejm/eejm1/fit.c', 'eejm/eejm1/common.c', 'eejm/eejm1/coord.c'],
   include_directories: [inc_np, include_directories('eejm/eejm1')],
   dependencies: [m_dep, gsl_dep],
   install: true,
   subdir: 'mpyricalspace'
)

py3.extension_module('eejm2',
   ['eejm/eejm2/eejm2.c', 'eejm/eejm2/eej_basis.c', 'eejm/eejm2/eej_calc.c',
    'eejm/eejm2/ndlinear.c', 'eejm/eejm2/common.c', 'eejm/eejm2/coord.c'],
   include_directories: [inc_np, include_directories('eejm/eejm2')],
   dependencies: [m_dep, gsl_dep],
   install: true,
   subdir: 'mpyricalspace'
)

py3.extension_module('eefm1',
   ['eefm1/eefm1.c', 'eefm1/eefm_basis.c', 'eefm1/eefm_calc.c', 'eefm1/ndlinear.c'],
   include_directories: [inc_np, include_directories('eefm1')],
   dependencies: [m_dep, gsl_dep],
   install: true,
   subdir: 'mpyricalspace'
)

# RTEEF -- the real-time (prompt-penetration) equatorial E-field model behind the
# Manoj & Maus PPEF web service, run locally from bundled ACE data (2001-2007).
py3.extension_module('ppeefm1',
   ['ppeefm1/ppeefm1.c', 'ppeefm1/rteef.c', 'ppeefm1/eefm_basis.c', 'ppeefm1/eefm_calc.c',
    'ppeefm1/ndlinear.c', 'ppeefm1/f107.c', 'ppeefm1/common.c', 'ppeefm1/lunar.c'],
   include_directories: [inc_np, include_directories('ppeefm1')],
   dependencies: [m_dep, gsl_dep],
   install: true,
   subdir: 'mpyricalspace'
)

# Keep the vendored IGRF Fortran source in sync with ngdc.noaa.gov before any
# igrfNN extension_module() below references it -- see igrf/sync_igrf.py's own
# docstring for what this does and does not touch, and why. Failures (no
# network, site unreachable) are swallowed inside the script itself, so this
# never blocks a build -- it just leaves the last-synced (committed) copy in
# place; check: false here is only a second line of defense.
sync_igrf_result = run_command(py3, files('igrf/sync_igrf.py'), check: false)
if sync_igrf_result.returncode() != 0
  warning('igrf/sync_igrf.py failed: ' + sync_igrf_result.stderr().strip())
endif

igrf14_module = custom_target('igrf14_module',
   output: 'igrf14fmodule.c',
   input: ['igrf/igrf14/igrf14f.pyf'],
   command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('igrf14f',
   ['igrf/igrf14/igrf14.f',
    igrf14_module,
   ],
   include_directories: [inc_np, inc_f2py],
   dependencies: fortranobject_dep,
   install: true,
   link_language: 'fortran',
   subdir: 'mpyricalspace',
)

igrf13_module = custom_target('igrf13_module',
   output: 'igrf13fmodule.c',
   input: ['igrf/igrf13/igrf13f.pyf'],
   command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('igrf13f',
   ['igrf/igrf13/igrf13.f',
    igrf13_module,
   ],
   include_directories: [inc_np, inc_f2py],
   dependencies: fortranobject_dep,
   install: true,
   link_language: 'fortran',
   subdir: 'mpyricalspace',
)

igrf12_module = custom_target('igrf12_module',
   output: 'igrf12fmodule.c',
   input: ['igrf/igrf12/igrf12f.pyf'],
   command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('igrf12f',
   ['igrf/igrf12/igrf12.f',
    igrf12_module,
   ],
   include_directories: [inc_np, inc_f2py],
   dependencies: fortranobject_dep,
   install: true,
   link_language: 'fortran',
   subdir: 'mpyricalspace'
)

igrf11_module = custom_target('igrf11_module',
   output: 'igrf11fmodule.c',
   input: ['igrf/igrf11/igrf11f.pyf'],
   command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('igrf11f',
   ['igrf/igrf11/igrf11.f',
    igrf11_module,
   ],
   include_directories: [inc_np, inc_f2py],
   dependencies: fortranobject_dep,
   install: true,
   link_language: 'fortran',
   subdir: 'mpyricalspace'
)

igrf10_module = custom_target('igrf10_module',
   output: 'igrf10fmodule.c',
   input: ['igrf/igrf10/igrf10f.pyf'],
   command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('igrf10f',
   ['igrf/igrf10/igrf10.f',
    igrf10_module,
   ],
   include_directories: [inc_np, inc_f2py],
   dependencies: fortranobject_dep,
   install: true,
   link_language: 'fortran',
   subdir: 'mpyricalspace'
)

igrf09_module = custom_target('igrf09_module',
   output: 'igrf09fmodule.c',
   input: ['igrf/igrf09/igrf09f.pyf'],
   command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('igrf09f',
   ['igrf/igrf09/igrf9.f',
    igrf09_module,
   ],
   include_directories: [inc_np, inc_f2py],
   dependencies: fortranobject_dep,
   install: true,
   link_language: 'fortran',
   subdir: 'mpyricalspace'
)

iri20_module = custom_target('iri20_module',
   output: 'iri20fmodule.c',
   input: ['iri/iri20/iri20f.pyf'],
   command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('iri20f',
   ['iri/iri20/irisub.for',
   'iri/iri20/irifun.for',
   'iri/iri20/cira.for',
   'iri/iri20/igrf.for',
   'iri/iri20/iridreg.for',
   'iri/iri20/iritec.for',
   'iri/iri20/iriflip.for',
   'iri/iri20/rocdrift.for',
   'iri/iri20/iri_batch.f90',
   iri20_module,
   ],
   include_directories: [inc_np, inc_f2py],
   dependencies: fortranobject_dep,
   install: true,
   link_language: 'fortran',
   subdir: 'mpyricalspace'
)

iri26_module = custom_target('iri26_module',
   output: 'iri26fmodule.c',
   input: ['iri/iri26/iri26f.pyf'],
   command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('iri26f',
   ['iri/iri26/irisub.for',
   'iri/iri26/irifun.for',
   'iri/iri26/cira.for',
   'iri/iri26/igrf.for',
   'iri/iri26/iridreg.for',
   'iri/iri26/iritec.for',
   'iri/iri26/iriflip.for',
   'iri/iri26/rocdrift.for',
   'iri/iri26/iri_batch.f90',
   iri26_module,
   ],
   include_directories: [inc_np, inc_f2py],
   dependencies: fortranobject_dep,
   install: true,
   link_language: 'fortran',
   subdir: 'mpyricalspace'
)

iri16_module = custom_target('iri16_module',
   output: 'iri16fmodule.c',
   input: ['iri/iri16/iri16f.pyf'],
   command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('iri16f',
   ['iri/iri16/irisub.for',
   'iri/iri16/irifun.for',
   'iri/iri16/cira.for',
   'iri/iri16/igrf.for',
   'iri/iri16/iridreg.for',
   'iri/iri16/iritec.for',
   'iri/iri16/iriflip.for',
   'iri/iri16/cosd_sind.for',
   'iri/iri16/iri_batch.f90',
   iri16_module,
   ],
   include_directories: [inc_np, inc_f2py],
   dependencies: fortranobject_dep,
   install: true,
   link_language: 'fortran',
   subdir: 'mpyricalspace'
)

iri12_module = custom_target('iri12_module',
   output: 'iri12fmodule.c',
   input: ['iri/iri12/iri12f.pyf'],
   command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('iri12f',
   ['iri/iri12/irisub.for',
   'iri/iri12/cira.for',
   'iri/iri12/igrf.for',
   'iri/iri12/iridreg.for',
   'iri/iri12/irifun.for',
   'iri/iri12/iritec.for',
   'iri/iri12/iriflip.for',
   'iri/iri12/iri_batch.f90',
   iri12_module,
   ],
   include_directories: [inc_np, inc_f2py],
   dependencies: fortranobject_dep,
   install: true,
   link_language: 'fortran',
   subdir: 'mpyricalspace'
)

# IRI-2007 and IRI-2001: no iriflip.for at these versions (introduced later).
iri07_module = custom_target('iri07_module',
   output: 'iri07fmodule.c',
   input: ['iri/iri07/iri07f.pyf'],
   command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('iri07f',
   ['iri/iri07/irisub.for',
   'iri/iri07/cira.for',
   'iri/iri07/igrf.for',
   'iri/iri07/iridreg.for',
   'iri/iri07/irifun.for',
   'iri/iri07/iritec.for',
   'iri/iri07/iri_batch.f90',
   iri07_module,
   ],
   include_directories: [inc_np, inc_f2py],
   dependencies: fortranobject_dep,
   install: true,
   link_language: 'fortran',
   subdir: 'mpyricalspace'
)

iri01_module = custom_target('iri01_module',
   output: 'iri01fmodule.c',
   input: ['iri/iri01/iri01f.pyf'],
   command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('iri01f',
   ['iri/iri01/irisub.for',
   'iri/iri01/cira.for',
   'iri/iri01/igrf.for',
   'iri/iri01/iridreg.for',
   'iri/iri01/irifun.for',
   'iri/iri01/iritec.for',
   'iri/iri01/iri_batch.f90',
   iri01_module,
   ],
   include_directories: [inc_np, inc_f2py],
   dependencies: fortranobject_dep,
   install: true,
   link_language: 'fortran',
   subdir: 'mpyricalspace'
)

#===============================================================================
# Coefficient / data files the models load at runtime -> mpyricalspace/_data/<model>/
# (without this the extension modules work only from a source checkout, not a wheel)
#===============================================================================
_pkg_data = py3.get_install_dir() / 'mpyricalspace' / '_data'

install_data(
  ['hwm/hwm14/hwm123114.bin', 'hwm/hwm14/dwm07b104i.dat', 'hwm/hwm14/gd2qd.dat'],
  install_dir: _pkg_data / 'hwm14',
)
install_data(
  ['hwm/hwm07/hwm071308e.dat', 'hwm/hwm07/dwm07b_104i.dat', 'hwm/hwm07/gd2qd.dat'],
  install_dir: _pkg_data / 'hwm07',
)
install_data(
  ['hltwim/hltwimnorth.dat', 'hltwim/hltwimsouth.dat', 'hltwim/gd2qd.dat'],
  install_dir: _pkg_data / 'hltwim',
)
install_data(
  ['eejm/eejm1/champ_mean_coeffs', 'eejm/eejm1/champ_stddev_coeffs'],
  install_dir: _pkg_data / 'eejm1',
)
install_data(
  ['eejm/eejm2/champ_mean_coeffs', 'eejm/eejm2/champ_stddev_coeffs',   # CHAMP (default)
   'eejm/eejm2/oer_mean_coeffs',   'eejm/eejm2/oer_stddev_coeffs',      # Orsted
   'eejm/eejm2/sac_mean_coeffs',   'eejm/eejm2/sac_stddev_coeffs'],     # SAC-C
  install_dir: _pkg_data / 'eejm2',
)
install_data(
  ['eefm1/champ_E_mean_coeffs', 'eefm1/champ_E_stddev_coeffs'],
  install_dir: _pkg_data / 'eefm1',
)
install_data(
  ['ppeefm1/TF.COF', 'ppeefm1/F107.txt',
   'ppeefm1/champ_E_mean_coeffs', 'ppeefm1/champ_E_stddev_coeffs',
   'ppeefm1/ace2001.dat.gz', 'ppeefm1/ace2002.dat.gz', 'ppeefm1/ace2003.dat.gz',
   'ppeefm1/ace2004.dat.gz', 'ppeefm1/ace2005.dat.gz', 'ppeefm1/ace2006.dat.gz',
   'ppeefm1/ace2007.dat.gz'],
  install_dir: _pkg_data / 'ppeefm1',
)

# apf107.dat/ig_rz.dat/CCIR/URSI come from iri/common/ (shared across every IRI
# version -- see iri/sync_iri.py's docstring) rather than each version's own folder.

# IRI-2016: IGRF-13 (igrf2020.dat), dgrf definitive to 2015
_iri16 = ['iri/common/apf107.dat', 'iri/common/ig_rz.dat', 'iri/iri16/igrf2020.dat', 'iri/iri16/igrf2020s.dat']
foreach m : range(11, 23)
  _iri16 += ['iri/common/ccir'  + m.to_string() + '.asc',
             'iri/common/ursi'  + m.to_string() + '.asc',
             'iri/iri16/mcsat' + m.to_string() + '.dat']
endforeach
foreach y : range(1945, 2016, 5)
  _iri16 += 'iri/iri16/dgrf' + y.to_string() + '.dat'
endforeach
install_data(_iri16, install_dir: _pkg_data / 'iri16')

# IRI-2020 (v2020.35) and IRI-2026: IGRF-14 (dgrf2020 definitive, igrf2025 provisional);
# IRI-2026 additionally ships ibp_emp_coeffs.dat for the Stolle et al. (2023) bubble model
foreach ver : ['iri20', 'iri26']
  _iri = ['iri/common/apf107.dat', 'iri/common/ig_rz.dat',
          'iri' / ver / 'igrf2025.dat', 'iri' / ver / 'igrf2025s.dat']
  foreach m : range(11, 23)
    _iri += ['iri/common/ccir'  + m.to_string() + '.asc',
             'iri/common/ursi'  + m.to_string() + '.asc',
             'iri' / ver / ('mcsat' + m.to_string() + '.dat')]
  endforeach
  foreach y : range(1945, 2021, 5)
    _iri += 'iri' / ver / ('dgrf' + y.to_string() + '.dat')
  endforeach
  if ver == 'iri26'
    _iri += 'iri/iri26/ibp_emp_coeffs.dat'
  endif
  install_data(_iri, install_dir: _pkg_data / ver)
endforeach

# IRI-2012: dgrf up to 2010, igrf2015 (not 2020), no mcsat storm model.
# apf107.dat/ig_rz.dat/ccir/ursi come from iri/common/ (shared across every
# IRI version -- see iri/sync_iri.py's docstring) rather than this version's
# own folder.
_iri12 = ['iri/common/apf107.dat', 'iri/common/ig_rz.dat']
foreach m : range(11, 23)
  _iri12 += ['iri/common/ccir' + m.to_string() + '.asc', 'iri/common/ursi' + m.to_string() + '.asc']
endforeach
_iri12 += ['iri/iri12/igrf2015.dat', 'iri/iri12/igrf2015s.dat']
foreach y : range(1945, 2011, 5)
  _iri12 += 'iri/iri12/dgrf' + y.to_string() + '.dat'
endforeach
install_data(_iri12, install_dir: _pkg_data / 'iri12')

# IRI-2007: dgrf up to 2005, igrf2010 -- same coefficient vintage as iri12's own
# (pre-sync) copy, just an older irisub.for/irifun.for (no iriflip.for at all).
_iri07 = ['iri/common/apf107.dat', 'iri/common/ig_rz.dat']
foreach m : range(11, 23)
  _iri07 += ['iri/common/ccir' + m.to_string() + '.asc', 'iri/common/ursi' + m.to_string() + '.asc']
endforeach
_iri07 += ['iri/iri07/igrf2010.dat', 'iri/iri07/igrf2010s.dat']
foreach y : range(1945, 2006, 5)
  _iri07 += 'iri/iri07/dgrf' + y.to_string() + '.dat'
endforeach
install_data(_iri07, install_dir: _pkg_data / 'iri07')

# IRI-2001: 2-digit dgrf/igrf file names (dgrf45.dat ... dgrf95.dat, dgrf00.dat,
# dgrf05.dat, igrf10.dat/igrf10s.dat) -- doesn't fit a clean numeric range, listed
# explicitly as published on the site.
_iri01 = ['iri/common/apf107.dat', 'iri/common/ig_rz.dat']
foreach m : range(11, 23)
  _iri01 += ['iri/common/ccir' + m.to_string() + '.asc', 'iri/common/ursi' + m.to_string() + '.asc']
endforeach
_iri01 += ['iri/iri01/igrf10.dat', 'iri/iri01/igrf10s.dat']
foreach y : ['45', '50', '55', '60', '65', '70', '75', '80', '85', '90', '95', '00', '05']
  _iri01 += 'iri/iri01/dgrf' + y + '.dat'
endforeach
install_data(_iri01, install_dir: _pkg_data / 'iri01')
