{-# LANGUAGE TemplateHaskell #-}
module Data.Chemistry.Wavefunction
( orb2AngMom
, angMom2Orb
, mixMOs
) where
import           Data.Chemistry.Types
import qualified Numeric.LinearAlgebra as BLAS

--------------------------------------------------------------------------------
-- Helpers for handling Wavefunction abstraction
--------------------------------------------------------------------------------
-- | Convert from orbital type to quantum number
orb2AngMom :: Char -> Int
orb2AngMom orb
    | orb == 's' || orb == 'S' = 0
    | orb == 'p' || orb == 'P' = 1
    | orb == 'd' || orb == 'D' = 2
    | orb == 'f' || orb == 'F' = 3
    | orb == 'g' || orb == 'G' = 4
    | orb == 'h' || orb == 'H' = 5
    | orb == 'i' || orb == 'I' = 6
    | orb == 'j' || orb == 'J' = 7
    | otherwise = error "not supported angular momentum"

-- | Angular momentum to orbital
angMom2Orb :: Int -> Char
angMom2Orb angmom
    | angmom == 0 = 's'
    | angmom == 1 = 'p'
    | angmom == 2 = 'd'
    | angmom == 3 = 'f'
    | angmom == 4 = 'g'
    | angmom == 5 = 'h'
    | angmom == 6 = 'i'
    | angmom == 7 = 'j'
    | otherwise = error "not supported angular momentum"


--------------------------------------------------------------------------------
-- calculations on Wavefunction
--------------------------------------------------------------------------------
-- give a list of weighted MOs and get back a new MO built of these and renormalized
mixMOs :: [(Double, MO)] -> MO
mixMOs lcl = newWMO
  where
    norm = sum . map fst $ lcl
    normCoeffs = map ((/ norm) . fst) lcl
    mos = map snd lcl :: [MO]
    rlcl = zip normCoeffs mos :: [(Double, MO)]
    nBF = length . BLAS.toList . head $ mos
    zeroMO = BLAS.fromList $ replicate nBF 0.0 :: MO
    newWMO = foldr (\(c, m) acc -> BLAS.fromList [c] * m + acc) zeroMO rlcl
