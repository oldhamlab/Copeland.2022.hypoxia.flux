# data.R

#' Cell count per DNA
#'
#' A data set containing the slopes to interpolate cell number from DNA
#' measurements.
#'
#' \describe{
#'   \item{cell_type}{
#'   `lf`    = lung fibroblasts \cr
#'   `pasmc` = pulmonary artery smooth muscle cells}
#'   \item{volume}{volume of lysis buffer used to extract DNA}
#'   \item{slope}{cells / ng DNA}
#'   }
"cells_per_dna"


#' Measurements for extracellular flux calculations
#'
#' A data set containing the interpolated metabolite concentrations, cell counts,
#' and evaporation volumes for extracellular flux determinations.
#'
#' \describe{
#'   \item{cell_type}{
#'   `lf` = lung fibroblasts \cr
#'   `pasmc` = pulmonary artery smooth muscle cells}
#'   \item{experiment}{
#'   `02` = 0.2% oxygen for hypoxia \cr
#'   `05` = 0.5% oxygen for hypoxia \cr
#'   `bay` = molidustat treatment \cr
#'   `05-bay` = 0.5% oxygen plus molidustat}
#'   \item{batch}{a group of biological replicates analyzed similarly}
#'   \item{date}{start date of an experiment}
#'   \item{metabolite}{name of measurement}
#'   \item{abbreviation}{abbreviated metabolite name}
#'   \item{detector}{
#'   `fld` = HPLC fluorescence detector \cr
#'   `mwd` = HPLC multi-wavelength detector \cr
#'   `enzyme` = enzymatic assay \cr
#'   `hplc` = OPD-derivatized HPLC detection \cr
#'   `lcms` = liquid chromatography-mass spectrometry assay \cr
#'   `picogreen` = fluorescence dye labeling of DNA}
#'   \item{type}{
#'   `cells` = conditioned medium \cr
#'   `empty` = unconditioned medium}
#'   \item{oxygen}{ambient oxygen level}
#'   \item{treatment}{
#'   `None` = no treatment \cr
#'   `DMSO` = 0.1% DMSO \cr
#'   `BAY` = 10 μM molidustat}
#'   \item{time}{hours}
#'   \item{well}{denotes technical replicates in each experiment}
#'   \item{conc}{measured in number for cells and μM for metabolites}
#'   \item{volume}{measured in mL, extrapolated from evaporation controls}
#'   \item{nmol}{metabolite mass}
#'   }
"flux_measurements"

#' Calculated fluxes
#'
#' A data set containing metabolite fluxes.
#'
#'  \describe{
#'   \item{metabolite}{name of measurement}
#'   \item{abbreviation}{abbreviated metabolite name}
#'   \item{cell_type}{
#'   `lf` = lung fibroblasts \cr
#'   `pasmc` = pulmonary artery smooth muscle cells}
#'   \item{experiment}{
#'   `02` = 0.2% oxygen for hypoxia \cr
#'   `05` = 0.5% oxygen for hypoxia \cr
#'   `bay` = molidustat treatment \cr
#'   `05-bay` = 0.5% oxygen plus molidustat}
#'   \item{batch}{a group of biological replicates analyzed similarly}
#'   \item{date}{start date of an experiment}
#'   \item{oxygen}{ambient oxygen level}
#'   \item{treatment}{
#'   `None` = no treatment \cr
#'   `DMSO` = 0.1% DMSO \cr
#'   `BAY` = 10 μM molidustat}
#'   \item{flux}{fmol / cell / h, positive fluxes indicate secretion, negative fluxes
#'   indicate uptake}
#'   }
"fluxes"

#' Growth rate parameters
#'
#' A data set containing the growth rate and X0 values from linear fitting of
#' cell count data.
#'
#' \describe{
#'   \item{cell_type}{
#'   `lf` = lung fibroblasts \cr
#'   `pasmc` = pulmonary artery smooth muscle cells}
#'   \item{experiment}{
#'   `02` = 0.2% oxygen for hypoxia \cr
#'   `05` = 0.5% oxygen for hypoxia \cr
#'   `bay` = molidustat treatment \cr
#'   `05-bay` = 0.5% oxygen plus molidustat}
#'   \item{batch}{a group of biological replicates analyzed similarly}
#'   \item{date}{start date of an experiment}
#'   \item{oxygen}{ambient oxygen level}
#'   \item{treatment}{
#'   `None` = no treatment \cr
#'   `DMSO` = 0.1% DMSO \cr
#'   `BAY` = 10 μM molidustat}
#'   \item{X0}{cell count at time 0}
#'   \item{mu}{growth rate per hour}
#'   }
"growth_rates"

#' Spontaneous degradation and accumulation rates
#'
#' A data set containing the rates of significant metabolite accumulation or
#' degradation in cell-free medium.
#'
#' \describe{
#'   \item{metabolite}{name}
#'   \item{oxygen}{ambient oxygen level}
#'   \item{treatment}{
#'   `None` = no treatment \cr
#'   `DMSO` = 0.1% DMSO \cr
#'   `BAY` = 10 μM molidustat}
#'   \item{k}{accumulation or decay rate per hour, positive values indicate decay}
#'   }
"k"

#' Quadrupole bias correction factors
#'
#' A data set containing experimentally determined correction factors to adjust
#' peak areas for isotopes detected in selected ion monitoring mode. Isotope peak
#' areas are multiplied by the correction factor. The corrected peak area is
#' subsequently used to determine the mass isotope distribution.
#'
#'  \describe{
#'    \item{batch}{correction factors associated with specific experiments}
#'    \item{metabolite}{name of measurement}
#'    \item{M}{isotope}
#'    \item{cf}{correction factor}
#'   }
"correction_factors"

#' Model reactions
#'
#' A data table of model reactions formatted for INCA analysis.
#'
#'  \describe{
#'   \item{index}{reaction identifier}
#'   \item{name}{reaction name}
#'   \item{equation}{equation for INCA input}
#'   \item{pathway}{metabolic pathway}
#'   }
"model_reactions"
