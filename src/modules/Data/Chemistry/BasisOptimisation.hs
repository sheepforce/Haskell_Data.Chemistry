module Data.Chemistry.BasisOptimisation
( getContrCoeff
, getContrCoeff_print
) where
--import qualified Numeric.LinearAlgebra as BLAS
--import qualified Data.Text as T
import           Data.Chemistry.BasisSet
import           Data.Chemistry.Molden
import           Data.Chemistry.Wavefunction
import           System.IO
import           Text.Printf
--import Data.List


-- get the contraction coefficients for a set of primitive GTOs (Basis functions)
-- give it:
--   a set of Basis functions, for which the contraction coefficient shall be determined "setOfBFs"
--   a set of molecular orbitals from which the contribution shall be calculated, e.g. three p orbitals
--     of a given shell, when calculating the contraction coefficients for p basis functions "setOfMOs"
--   if the resulting contraction coefficients should be renormalized or taken plain from the calculation "renorm"
--   the molden file, from which the calculation is going to be done
getContrCoeff :: [Int] -> [Int] -> Bool -> Molden -> [Double]
getContrCoeff setOfBFs setOfMOs renorm molden
    | all (== ((head angMomOfSelBFs) :: Int)) (tail angMomOfSelBFs) == True && renorm == False = contrCoeff
    | all (== ((head angMomOfSelBFs) :: Int)) (tail angMomOfSelBFs) == True && renorm == True = contrCoeff_renorm
    | otherwise = error "trying to combine basis functions of different angular momentum"
    where
        -- for making sure no nonsensical combinations of basis functions were selected
        angMomOfAllBFs = concat . map (map basfun_angular) . basfuns $ molden
        angMomOfSelBFs = [angMomOfAllBFs !! i | i <- setOfBFs]

        -- gives [[Double]], where outer layer is the MO, from which coefficients were determined
        -- and inner layer is the basis function
        moCoeffsByMOandBF = [[getMOcoeffsForBFinOrb j i (head angMomOfSelBFs) molden | i <- setOfBFs] | j <- setOfMOs]

        -- now properly average them together (over the MOs)
        contrCoeff = [(sum $ map (!! i) moCoeffsByMOandBF) / (fromIntegral $ length moCoeffsByMOandBF) | i <- [0 .. ((length $ head moCoeffsByMOandBF) - 1)]]

        -- renormalize contraction coefficient
        renormN = 1.0 / (sum contrCoeff)
        contrCoeff_renorm = map (* renormN) contrCoeff

getContrCoeff_print :: Handle -> [Int] -> [Int] -> Bool -> Molden -> IO()
getContrCoeff_print handle setOfBFs setOfMOs renorm molden = do
    -- calculate the contraction coefficients
    let contrCoeff = getContrCoeff setOfBFs setOfMOs renorm molden
        angMomOfAllBFs = concat . map (map basfun_angular) . basfuns $ molden
        angMomOfSelBFs = [angMomOfAllBFs !! i | i <- setOfBFs]
        angMomOfSelBF = head angMomOfSelBFs


    -- print them
    putStrLn $ "contraction coefficients for"
    putStrLn $ "  angular momentum : " ++ show (angMom2Orb angMomOfSelBF)
    putStrLn $ "  basis functions  : " ++ show setOfBFs
    putStrLn $ "  from the MOs     : " ++ show setOfMOs
    putStrLn $ "  renormalized     : " ++ show renorm
    putStrLn $ ""
    mapM_ (hPrintf handle "%9.7f\n") contrCoeff
