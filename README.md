# PUMA
The Positional Update and Matching Algorithm (PUMA) arose from the need for accurate foreground models for calibration of low radio frequency interferometric data. It was conceived to leverage the positional accuracies of higher frequency catalogues, whilst using lower frequency catalogues to obtain accurate spectral energy distributions (SEDs). The difficulty in cross matching multiple surveys is the varying instrument parameters for each survey such as angular resolution and sensitivity, as well changing morphologies of sources with frequency. PUMA initially uses a bayesian positional probability match based on (Tamas and Budavari, 2008), to identify any possible matches. It then uses the fact that at low radio frequencies, most SEDs are expected to be appoximately linear in log log space. PUMA uses this when higher resolution catalogues report multiple possible matched components, by adding the flux of the higher resolution sources and testing for linearity.

### 30/06/2015
PUMA is currently undergoing a bit of a face lift, and as such catalogue_algorithm.pdf is somewhat out of date. Hopefully I'll be able to update it soon. The algorithms are the same, but the example scripts within are partially obsolete. The scripts within the dir examples are up to date however, so you should be able to run them after following the steps in install_notes.txt. If you have any questions, feel free to email me at jline@student.unimelb.edu.au.

Thanks!

Jack Line
