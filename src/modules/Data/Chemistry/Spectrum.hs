module Data.Chemistry.Spectrum
( gauss
, convolutionSum
, oscStrength2Epsilon
, eV2nm
, nm2eV
) where

-- | A gaussian line shape as defined in Gaussian whitepaper http://gaussian.com/uvvisplot/ .
-- | The definition used here differs from the Wikipedia version, in the prefactor. Here the area
-- | is not normalised but the height, so that every peak can be safely multiplied with the
-- | oscillator strength.
-- |   σ -> Standard deviation for the line width. This can be converted to FWHM by 2*sqrt(ln 2),
-- |        instead of the 2*sqrt(2 * ln 2), coming from the σ instead of σ^2, i think
-- |   fwhm -> Full width at Half Maximum, directly converted to standard deviation
-- |   eEval -> Energy in electron Volt at which to obtain the value of the broadened spectrum
-- |   ePeak -> Energy in electron Volt for the peak position (x-Coordinate)
-- |   fOsc -> Dimensionaless oscillator strength as in the QC output file
gauss :: Floating a => a -> (a, a) -> a -> a
gauss fwhm (ePeak, fOsc) eEval = fOsc * exp(- ((eEval - ePeak) / σ)**2)
  where
    σ = fwhm / (2 * sqrt(log 2))

-- | Broaden a complete stick spectrum by a given (parametrised) broadening function.
-- |   convFunc -> Function taking pairs of (peak position in electron Volt, oscillator strength)
-- |               and a grid point, on which to evaluate the value
-- |   peaks -> A list of pairs of (peak position in electron Volt, oscillator strength)
-- |   grid -> A list of grid points (x-Values in electron Volt)
-- | Result will be a list of (grid points in electron Volt, broadened oscillator strength at pos.)
convolutionSum ::  Floating a => ((a, a)-> a -> a) -> [(a, a)] -> [a] -> [(a, a)]
convolutionSum convFunc peaks grid = gridResult
  where
    convFunc' grid' peaks' = convFunc peaks' grid'
    gridPointVal p = sum . map (convFunc' p) $ peaks
    gridVal = map gridPointVal grid
    gridResult = zip grid gridVal

-- | Convert from oscillator strength (dimensionless) to absorption coefficients (litre / (mol cm))
-- | The prefactor is taken from the Gaussian whitepaper, and the nominator multiplicator comes
-- | from the prefactor expecting σ to be in cm^-1, but i give it in electron Volts instead
oscStrength2Epsilon :: Floating a => a -> a -> a
oscStrength2Epsilon fwhm fOsc = 1.30062974e8 * (fOsc) / (σ * 8065.54400545911)
  where
    σ = fwhm / (2 * sqrt(log 2))

-- | Convert from electron Volt to nano metre
eV2nm :: Floating a => a -> a
eV2nm x = 1239.84197386209 / x

-- | Convert from nano metre to electron Volt
nm2eV :: Floating a => a -> a
nm2eV = eV2nm
