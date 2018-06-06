{-# LANGUAGE TemplateHaskell #-}
module Data.Chemistry.Molden
( n_e
, n_BF
, angMom_BF
, bf2aos_sph
, basis2AtomicOrbitals
, moOrderingInAngMom
, moOrderingInAtoms
, moOrderingInBF
, getMOcoeffsForBFinOrb
, mixMMOs
) where
import           Data.Chemistry.Types
import           Data.Chemistry.Wavefunction
import           Data.List
import           Lens.Micro.Platform
import qualified Numeric.LinearAlgebra       as BLAS


--------------------------------------------------------------------------------
-- working with molden files
--------------------------------------------------------------------------------
-- | Get the overall number of electrons in the molecule
n_e :: Molden -> Double
n_e molden = sum . map _moldenMO_occup $ (molden ^. molden_mos)

-- | Get the overall number of basis functions in a molden file
n_BF :: Molden -> Int
n_BF molden = length . concat $ molden ^. molden_basfuns

-- | Atomwise list of angular momentums for the basis functions
angMom_BF:: Molden -> [[Int]]
angMom_BF molden = map (map _basfun_angular) $ _molden_basfuns molden

-- | Give in a basis function and you will get atomic orbitals from it in molden ordering
bf2aos_sph :: BasFun -> [AO]
bf2aos_sph bf
    -- this is the molden spherical ordering for s
    | bf_angmom == 0 = [ AO { _ao_n = 0
                            , _ao_l = bf_angmom
                            , _ao_m = i
                            , _ao_r = cgto
                            } | i <- [0]
                       ]
    -- this is the molden spherical ordering for p
    | bf_angmom == 1 = [ AO { _ao_n = 0
                            , _ao_l = bf_angmom
                            , _ao_m = i
                            , _ao_r = cgto
                            } | i <- [(-1), 0 , 1]
                       ]
    -- this is the molden spherical ordering for d
    | bf_angmom == 2 = [ AO { _ao_n = 0
                            , _ao_l = bf_angmom
                            , _ao_m = i
                            , _ao_r = cgto
                            } | i <- [0, 1, (-1), 2, (-2)]
                       ]
    -- this is the molden spherical ordering for f
    | bf_angmom == 3 = [ AO { _ao_n = 0
                            , _ao_l = bf_angmom
                            , _ao_m = i
                            , _ao_r = cgto
                            } | i <- [0, 1, (-1), 2, (-2), 3, (-3)]
                       ]
    -- this is the molden spherical ordering for f
    | bf_angmom == 4 = [ AO { _ao_n = 0
                            , _ao_l = bf_angmom
                            , _ao_m = i
                            , _ao_r = cgto
                            } | i <- [0, 1, (-1), 2, (-2), 3, (-3), 4, (-4)]
                       ]
    | bf_angmom == 5 = [ AO { _ao_n = 0
                            , _ao_l = bf_angmom
                            , _ao_m = i
                            , _ao_r = cgto
                            } | i <- [0, 1, (-1), 2, (-2), 3, (-3), 4, (-4), 5, (-5)]
                       ]
    | otherwise = error "Molden does not support functions higher than g"
    where bf_angmom = _basfun_angular bf
          cgto = _basfun_radial bf

-- | from the basis functions expand in the AOs
-- | ordering of the list
-- |   atom
-- |     Hauptquantenzahl, Nebequantenzahl (e.g. 1s / 2s / 3d / 4f / 5f)
-- |       Magnetquantenzahl
basis2AtomicOrbitals :: Molden -> [[[AO]]]
basis2AtomicOrbitals molden = mol_ao
    where
        -- the basis functions atomwise
        mol_bf = molden ^. molden_basfuns

        -- get angular momentums of the basis functions atomwise
        mol_ao = map (map bf2aos_sph) mol_bf

-- | The order of the MOs in the basis of the AOs (giving angular momentum here)
moOrderingInAngMom :: Molden -> [Int]
moOrderingInAngMom molden = concat . concat . map (map (map _ao_l)) . basis2AtomicOrbitals $ molden

-- | The order of the MOs in the basis of the AOs (giving in corresponding atom)
moOrderingInAtoms :: Molden -> [Int]
moOrderingInAtoms molden =  fillLLwithPos . map concat . basis2AtomicOrbitals $ molden

-- | The order of the MOs in the basis of the AOs (given in corresponding basis functions)
moOrderingInBF :: Molden -> [Int]
moOrderingInBF molden = fillLLwithPos . concat . map (map (map _ao_l)) . basis2AtomicOrbitals $ molden

-- | Get MO coefficients (sum over AOs with same n and l but different m)
-- | for a specific basis function in a molecular orbital
getMOcoeffsForBFinOrb :: Int -> Int -> Int -> Molden -> Double
getMOcoeffsForBFinOrb mo_number bf_number angMomOfBF molden
    |  listIndexForBF `isInfixOf` listIndexForAngMom = sum coeffs_for_BF
    | otherwise = error "the desired angular momentum does not belong to the specified basis function"
    where
        mocoeffs = _moldenMO_coeffs $ (molden ^. molden_mos) !! mo_number
        bfIndexList = moOrderingInBF molden
        angMomIndexList = moOrderingInAngMom molden
        listIndexForBF = findIndices (== bf_number) bfIndexList
        listIndexForAngMom = findIndices (== angMomOfBF) angMomIndexList
        coeffs_for_BF = [mocoeffs BLAS.! i | i <- listIndexForBF]


--------------------------------------------------------------------------------
-- Manipulation of Molden files
--------------------------------------------------------------------------------
-- | arbitrary linear combinations of MOs from Molden files, that produces a new
-- | Molden file with a single orbital, which is the new linearily combined one
mixMMOs :: Molden -> [(Double, Int)] -> Molden
mixMMOs molden wMOList = newMolden
  where
    orbWeights = map fst wMOList
    orbIndices = map snd wMOList
    orbs2Mix = [ map _moldenMO_coeffs (molden ^. molden_mos) !! i | i <- orbIndices ]
    linearCombList = zip orbWeights orbs2Mix
    newMOCoeffs = mixMOs linearCombList
    newMMO = MoldenMO
      { _moldenMO_sym    = "c1"
      , _moldenMO_energy = 0.0
      , _moldenMO_spin   = Alpha
      , _moldenMO_occup  = 2.0
      , _moldenMO_coeffs = newMOCoeffs
      }
    newMolden = molden & molden_mos .~ [ newMMO ]


--------------------------------------------------------------------------------
-- non specific helper functions
--------------------------------------------------------------------------------
-- | take a nested list and fills the inner entries with its position in the outer index, than flatten
-- | [[1,2,3,4], [5,6,7]] -> [0,0,0,0,1,1,1]
fillLLwithPos :: [[a]] -> [Int]
fillLLwithPos a = concat [Prelude.take (length $ a !! i) $ repeat i | i <- [0 .. ((length a) - 1)]]
