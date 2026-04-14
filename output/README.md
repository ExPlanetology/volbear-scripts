## Data Repository & Output Structure

**NOTE:** The GitHub repository ((https://github.com/ExPlanetology/volbear-scripts)) contains the processing scripts only. The full atmospheric model outputs are hosted on **Zenodo**: 10.5281/zenodo.19571079

---

The simulation data is organized into the following directories, corresponding to the grids and figures presented in the associated manuscript:

### 1. Primary Atmosphere Grids
* **`grid/`**
  Contains fully equilibrated atmospheres for the primary model grid. These simulations are restricted to a fully molten **55 Cnc e** analogue (8 Earth masses) across the following parameter space:
  * **Psi (Ψ)**: [-1, 0, 1]
  * **Xi (ξ)**: [0, 0.1, 0.2, 0.5, 1]
  * **ΔIW**: [-6, -3, 0, 3, 6]
  * **T_irr**: [1000, 1500, 2000, 2500] K
* **`grid_with_aCeTY/`**
  A subset of the primary grid recalculated with C2H2 opacities included. (Reference: Appendix J, Figs. J.1 & J.2).
* **`grid_with_CIA/`**
  Selected cases incorporating Collision Induced Absorption (CIA). (Reference: Fig. B.1).

### 2. Sensitivity & Series Studies
* **`mass_series/`**
  Atmospheric profiles modeled as a function of planetary mass. Parameters are restricted to Ψ = 0 with specific variations in ΔIW [-3, 3], ξ [0, 1], and T_irr [1250, 2500] K. (Used in Fig. 9).
* **`tint_series/`**
  Simulations varying the intrinsic temperature (T_int) for 8 Earth-mass planets, following a similar parameter restriction as the mass series. (Reference: Appendix F).

### 3. Specialized Models
* **`MAI_only/`**
  Model outputs for "Magma-Atmosphere Interaction only" cases.
  * **`batch/`**: Batch calculations utilized in Figs. G.1 and G.2.
  * **`fugacity_series/`**: Outputs specifically generating the results in Figs. 2 and 3.
