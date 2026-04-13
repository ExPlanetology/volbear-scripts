# volbear

Scripts and models used in "Volatile-bearing mineral atmospheres of hot rocky exoplanets as probes of interior state and composition".

## Description

This repository archives the scripts and models used in the study "Volatile-bearing mineral atmospheres of hot rocky exoplanets as probes of interior state and composition." It is intended as documentation rather than standalone software.

The study uses a modified version of the MAGMA code (Fegley & Cameron 1987; Fegley & Schaefer 2009), herein referred to as muspell (see Seidler et al. 2024 for details). MAGMA predates modern permissive software licensing, and its current licensing status — including that of all derivatives such as muspell — is unclear; consequently muspell is not included in this repository.

## License

This project is licensed under the Affero General Public License (AGPL). See the [LICENSE](LICENSE) file for more details.

### Figures

The figures included in this repository are licensed under a Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0).

* You are free to share and adapt the images, but not for commercial purposes.
* You must give appropriate credit, provide a link to the license, and indicate if changes were made.
* You may not use the material for commercial purposes, including but not limited to use in any product, service, or promotion intended for sale or profit.
* For more details on the terms of use, please refer to the full CC BY-NC 4.0 license.

## Installation

The project is installable by standard means, e.g. pip, uv, poetry, etc. Note, however, that it cannot be executed without an installation of `muspell`, which is not publicly available (see **Description** above).

Opacity files are required. These can be build manually (see https://github.com/ExPlanetology/phaethon).

An environment variable (OPAC_PATH) has to point to the path containing the opacities: 
```
    export OPAC_PATH=/path/to/opacities/
```

Further, it is crucial for a CUDA-capable GPU to be present, and for the code to be compiled for the correct architecture, which has to be done manually (see https://github.com/ExPlanetology/phaethon/blob/main/README.md)
