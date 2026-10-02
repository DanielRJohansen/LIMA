# LIMA
LIMA is a GPU accelerated molecular dynamics engine, with tools to build, solvate, energy-minimize, simulate and
render molecular systems.

## Install
Download the latest release for Windows, Ubuntu/Debian, Arch or other Linux from the
[Releases page](https://github.com/DanielRJohansen/LIMA/releases), and run `lima --help` to get started.

LIMA requires an NVIDIA GPU of the RTX 40-series or newer (or H100, B200 and similar), with a recent NVIDIA driver.

## License
LIMA is free for small companies, for noncommercial use and for academia, and can be evaluated by anyone.
It is available under your choice of three licenses, see [LICENSE.txt](LICENSE.txt):

- **PolyForm Small Business 1.0.0**: any use by companies with fewer than 100 people and less than 1,000,000 USD
  revenue.
- **PolyForm Noncommercial 1.0.0**: noncommercial use, and use by universities, public research organizations,
  charities and government institutions.
- **PolyForm Free Trial 1.0.0**: evaluation by anyone for less than 32 consecutive days.

Other use, such as by a larger company beyond evaluation, requires a commercial license.
Contact daniel@lima-dynamics.com.

LIMA includes third-party software under its own licenses, see [THIRD_PARTY_NOTICES.txt](THIRD_PARTY_NOTICES.txt).

## Building from source
- Windows: Visual Studio 2022 with the C++ workload and the CUDA toolkit. Open the folder in Visual Studio, or
  configure with CMake and Ninja.
- Linux: run `install.sh`, which builds LIMA and installs it to /usr/bin and /usr/share/LIMA.

Releases are built with `distribution/release.bat`, see [distribution/README.md](distribution/README.md).


## LIMA would not be possible without scientific contributions of many researchers. Citations for resources below:

### LIPIDS
Lipid files supplied by Stockholm Lipids:
http://www.fos.su.se/~sasha/SLipids/Downloads.html
Dr. Joakim Jämbeck
Prof. Alexander Lyubartsev


Saturated PC lipids:
Derivation and Systematic Validation of a Refined All-Atom Force Field for Phosphatidylcholine Lipids,            
Joakim P. M. J�mbeck, Alexander P. Lyubartsev, J. Phys. Chem. B, 2012, 116 (10), 3164-3179
DOI: 10.1021/jp212503e


POPC, DOPC, SOPC, DOPE, POPE and similar:
An Extension and Further Validation of an All-Atomistic Force Field for Biological Membranes, Joakim P. M. J�mbeck, Alexander P. Lyubartsev, J. Chem. Theory Comput., 8 (8), 2938-2948 DOI: 10.1021/ct300342n

PS, PG, SM lipids and Cholesterol:
Another Piece of the Membrane Puzzle: Extending Slipids Further,
Joakim P. M. J�mbeck, Alexander P. Lyubartsev, J. Chem. Theory Comput., 9 (1), 774-784 (2013) 
DOI: 10.1021/ct300777p


Polyunsaturated lipids:
Extension of the Slipids Force Field for Polyunsaturated Lipids,
Inna Ermilova, Alexander P. Lyubartsev, J. Phys. Chem. B, 120 (50), pp 12826�12842 (2016) 
DOI: 10.1021/acs.jpcb.6b05422


Slipids update 2020:
Optimization of Slipids Force Field Parameters Describing Headgroups of Phospholipids,
Fredrik Grote, Alexander P. Lyubartsev, J. Phys. Chem. B, 124 (50), pp 8784-8793 (2020) 
DOI: 10.1021/acs.jpcb.0c06386 









