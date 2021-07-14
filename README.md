# PUMA

The Positional Update and Matching Algorithm (PUMA) arose from the need for accurate foreground models for calibration of low radio frequency interferometric data. It was conceived to leverage the positional accuracy of higher frequency catalogues, whilst using lower frequency catalogues to obtain accurate spectral energy distributions (SEDs). The difficulty in cross matching multiple surveys is the varying instrument parameters for each survey such as angular resolution and sensitivity, as well changing morphologies of sources with frequency. PUMA initially uses a Bayesian positional probability match based on [Budavari & Szalay, 2008](https://iopscience.iop.org/article/10.1086/587156), to identify any possible matches. It then uses the fact that at low radio frequencies, most SEDs are expected to be approximately linear in log log space. PUMA uses this when higher resolution catalogues report multiple possible matched components, by adding the flux of the higher resolution sources and testing for linearity.

For details on the algorithms and an example application, see [Line et al. 2017](https://doi.org/10.1017/pasa.2016.58).

## Dependencies

### stilts

Installation/download instructions [for `stilts` here](http://www.star.bristol.ac.uk/~mbt/stilts/)

I've removed the bundled `topcat` and `stilts` versions that used to come with `PUMA` as it was quite rightly pointed out they have aged with time, and could interfere with system versions. And `topcat` isn't even used by `PUMA` (what was I thinking?) So now you need to install `stilts` yourself. If you use something like Ubuntu, you'll need to install `java` via something like

```
sudo apt install default-jre
```

I installed `stilts` via the jar file and bash script, and put it in my `/usr/local/bin` so you can find it on the command line. To do that I did:

```bash
cd /usr/local/bin
sudo wget http://www.star.bristol.ac.uk/~mbt/stilts/stilts.jar .
sudo wget http://www.star.bristol.ac.uk/~mbt/stilts/stilts
chmod +x stilts

```

### python
There are a few `python` dependencies, which are listed in `requirements.txt`, which you can install with something via

```
pip --install -r requirements.txt
```

## `PUMA` Installation

To download, simply git clone with the following command:

```
git clone https://github.com/JLBLine/PUMA.git
```

Once downloaded, add these two lines to your `~/.bashrc` to find PUMA scipts:

```
export PUMA_DIR=*path-to-you-PUMA-directory*
export PATH=$PUMA_DIR/scripts:$PATH
```

## Example run

The example bash scripts in `$PUMA_DIR/examples` use the `$PUMA_DIR` variable, so they will work out of the box if you follow the installation notes above, as well as heading to `$PUMA_DIR/original_cats`, and unpacking the catalogues:

```
tar -vxjf vlssr_names.fits.tar.bz2
tar -vxjf vizier_mwacs.fits.tar.bz2
tar -vxjf vizier_mrc.fits.tar.bz2
tar -vxjf sumss_names.fits.tar.bz2
tar -vxjf vizier_nvss.fits.tar.bz2
```

Once this is done, typing
```
source timing.sh
```
at the terminal will run `PUMA` to match MWACS to VLSSr, MRC, SUMSS and NVSS.

After this, running

```
source plot_outcomes.sh
```

in the same directory will plot some of the outcomes. If you run `plot_outcomes.py --help`, the options will be summarised (there are quite a few!).

Thanks!

Jack Line
