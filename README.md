# spm12-dartel

Code for preprocessing of functional and structural MRI data into standardized MNI space using SPM12 and DARTEL. 

* Can be used with only one structural scan (e.g. either T1 MPRAGE or T2 matched-bandwidth)
* Can be used with two structural scans (e.g. T1 MPRAGE *and* T2 matched-bandwidth). Secondary scan (e.g. MBW) used as intermediary for coregistering functionals to primary structural (e.g. MPRAGE)

<b>Instructions:</b>

Call only the <b>wrapper</b> script as it will call the <b>run</b> function in a parfor loop. All user-editable parameters are in the epynomous section of the wrapper. Other sections of the wrapper script and run function shouldn't be edited unless you know what you're doing.

A "runStatus" struct containg each subject's pre-dartel status will be saved in the folder specified in "batchDir". The matlab workspace after pre-dartel will also be saved in "batchDir", you can use this to re-run DARTEL without re-running pre-dartel. A text log of the matlab console output will be saved for predartel & dartel in the "batchDir" folder. All pre-dartel and DARTEL matlabbatches will be saved in "batchDir"

Final image ouputs for further analyses:
* rBOLD[run name].nii = resliced, native-space BOLD images (e.g. native-space single-subject analyeses)
* wrBOLD[run name].nii = resliced, MNI-space & unsmoothed BOLD images (e.g. connectivity or MVPA analyses)
* swrBOLD[run name].nii = resliced, MNI-space & smoothed BOLD images (e.g. normal group-level analyeses)
* rp_BOLD[run name].nii = motion parameters to include in level1 model as nuisance regressors
* wm[MPRAGE/structural name].nii = bias-corrected MNI-space structural image
* wc1[MPRAGE/structural name].nii = MNI-space grey matter segmentation image (can be used as mask/ROI; signal from which typically used for functional connectivity analyses, etc)
* wc2[MPRAGE/structural name].nii = MNI-space white matter segmentation image (can be used as mask/ROI; signal from which typically regressed out during functional connectivity analyses, etc)
* wc3[MPRAGE/structural name].nii = MNI-space cerebrospinal fluid (also eyes/sinuses maybe) segmentation image (can be used as mask/ROI; signal from which typically regressed out during functional connectivity analyses, etc)

<b>Algorithm when using only one structural scan:</b>
1) Realign & reslice functionals to mean functional (pre-dartel; parfor parallelization)
2) Co-register structural to mean functional (pre-dartel; parfor parallelization)
3) Segment & bias-correct structural, generate segment params for DARTEL  (pre-dartel; implicit multithreading from here)
4) Create DARTEL templates & generate deformation fields for MNI normalization
5) Normalize & smooth functionals to MNI space via DARTEL (implicit parallelization from here)
6) Normalize bias-corrected structural to MNI space via DARTEL
7) Normalize grey matter (C1) segmentation to MNI space via DARTEL
8) Normalize white matter (C2) segmentation to MNI space via DARTEL
9) Normalize CSF (C3) segmentation to MNI space via DARTEL

<b>Algorithm when using MPRAGE & MBW:</b>
1) Realign & reslice functionals to mean functional (pre-dartel; parfor parallelization)
2) Co-register MBW to mean functional (pre-dartel; parfor parallelization)
3) Co-register MPRAGE to MBW (pre-dartel; parfor parallelization)
4) Segment & bias-correct structural, generate segment params for DARTEL  (pre-dartel; implicit multithreading from here)
5) Create DARTEL templates & generate deformation fields for MNI normalization
6) Normalize & smooth functionals to MNI space via DARTEL (implicit parallelization from here)
7) Normalize bias-corrected MPRAGE to MNI space via DARTEL
8) Normalize grey matter (C1) segmentation to MNI space via DARTEL
9) Normalize white matter (C2) segmentation to MNI space via DARTEL
10) Normalize CSF (C3) segmentation to MNI space via DARTEL
