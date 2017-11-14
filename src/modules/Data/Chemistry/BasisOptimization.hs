import qualified Data.Text as T
import Data.Chemistry.Molden
import Data.Chemistry.Wavefunction
import Data.List
import qualified Numeric.LinearAlgebra as BLAS


-- get MO coefficients for a specific basis function in a molecular orbital
getMOcoeffsForBFinOrb :: Int -> Int -> Int -> Molden -> [Double]
getMOcoeffsForBFinOrb mo_number bf_number angMomOfBF molden
    | isInfixOf listIndexForBF listIndexForAngMom == True = coeffs_for_BF
    | otherwise = error "the desired angular momentum does not belong to the specified basis function"
    where
        mocoeffs = coeffs $ (mos molden) !! mo_number
        bfIndexList = moOrderingInBF molden
        angMomIndexList = moOrderingInAngMom molden
        listIndexForBF = findIndices (== bf_number) bfIndexList
        listIndexForAngMom = findIndices (== angMomOfBF) bfIndexList
        coeffs_for_BF = [mocoeffs BLAS.! i | i <- listIndexForBF]
        
