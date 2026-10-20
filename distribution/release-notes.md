## What's new

### Faster
- **Molecular dynamics is about 1.8x faster** on large systems (STMV benchmark, ~1M atoms: ~7.6 -> ~4.1 ms/step on an RTX 40-series), from reworked nonbonded, bonded, PME, integration and neighbor-search kernels. Small systems gain too.
- **Energy minimization** uses the same faster nonbonded kernel.
- **Faster startup**: a small system goes from command to first step in ~45 ms instead of ~115 ms, and batches of simulations start sooner.

### More accurate
- **PME is now standard smooth PME**: long-range electrostatic force error dropped from ~17% to ~0.3%. Results with electrostatics will differ from 1.3.
- **The thermostat runs on its own schedule.** It used to depend on trajectory logging, so with logging off it never ran.
- **Solvation keeps water 0.35 nm away from the solute**, so membranes no longer get water in their hydrophobic core.

### New in the CLI
- **GPU commands share the GPU**: `lima mdrun`, `em` and `buildmembrane` started with `--async` run on a background `lima` that batches them on the GPU; `lima wait` collects the results. Without `--async` they run in-process. Use `lima daemon status|stop` to manage it.
- **`--async`** returns immediately, and `lima wait` collects the results. With `-d`, a batch is shown tiled in one window.
- **`lima mdrun --set KEY=VALUE ...`** overrides parameters of the sim_params file.
- **`initial_temperature`** in sim_params draws starting velocities.
- **`lima em`** reports its steps, final max force and batch size.
- **Windows**: on first run, lima offers to add itself to PATH and to open `.gro` files. `lima setup` does this later, and `lima setup --undo` reverts it.

### Fixes
- Inputs that exceed the engine's internal capacities now stop with an error instead of silently corrupting the simulation.
- Invalid input is rejected with a clear error: NaN positions, invalid dt, cutoff or masses, and overlapping particles in MD.
- Energy minimization no longer reports convergence on NaN forces, and handles straight angles and exactly overlapping particles.

## Downloads

| Platform | File |
|---|---|
| Windows 10/11 (x64) | `lima-@VERSION@-windows-x64.zip`: unzip anywhere and run `lima`, which offers to add itself to PATH |
| Ubuntu / Debian | `lima_@VERSION@_amd64.deb`: `sudo apt install ./lima_@VERSION@_amd64.deb` |
| Arch | `PKGBUILD`: `makepkg -si` in a folder containing it |
| Other Linux (x86_64) | `lima-@VERSION@-linux-x86_64.tar.gz`: extract anywhere and run `bin/lima` (needs glfw, OpenGL and oneTBB) |

## Requirements

An NVIDIA GPU of the RTX 40-series or newer (or H100, B200 and similar), with a recent NVIDIA driver (CUDA 13 capable).

Run `lima --help` to get started.

## License

Free for small companies (fewer than 100 people and under 1M USD revenue), for noncommercial use and for academia, and free to evaluate for 31 days for anyone, under the PolyForm Small Business, Noncommercial and Free Trial licenses. Other use requires a commercial license: daniel@lima-dynamics.com. See `LICENSE.txt` in the download.
