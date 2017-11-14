module Data.Chemistry.Molden
( MoldenCoord(..)
, MMO(..)
, Molden(..)
, moldenParser
, n_e
, n_BF
, angMom_BF
, bf2aos_sph
, basis2AtomicOrbitals
, moOrderingInAngMom
, moOrderingInAtoms
, moOrderingInBF
, getMOcoeffsForBFinOrb
) where
import qualified Numeric.LinearAlgebra as BLAS
import qualified Data.Text as T
import Control.Applicative
import Data.Attoparsec.Text.Lazy
--import Data.Chemistry.XYZ
import Data.Chemistry.Wavefunction
import Data.Chemistry.BasisSet
import Data.List
import Data.Maybe()


{- ############################ -}
{- define data types for Molden -}
{- ############################ -}

-- units for geometry may either be Angstrom or Bohr
data Units  = Angstrom | Bohr deriving Show

-- a data type for the geometry informations in a molden [Atoms] block
data MoldenCoord = MoldenCoord { moco_element :: String
                               , moco_natom :: Int
                               , moco_nproton :: Int
                               , moco_coord :: (Double, Double, Double)
                               } deriving Show

-- extended data type for MOs in a Molden-File
data MMO = MMO { sym :: String
               , energy :: Double
               , spin :: Spin
               , occup :: Double
               , coeffs :: MO
               } deriving Show

-- contents of a Molden file
data Molden = Molden { -- atoms contains the geometry
                       atoms :: (Units, [MoldenCoord])
                       -- the basis functions, outer for the atom, inner for gto
                     , basfuns :: [[BasFun]]
                     , mos :: [MMO]
                     } deriving Show


{- ###### -}
{- Parser -}
{- ###### -}

moldenParser :: Parser Molden
moldenParser = do
    -- it all began with "[Molden Format]" ...
    _ <- string $ T.pack "[Molden Format]"    
    
    -- parse the "[Atoms]" block, containing the geometry informations
    skipSpace
    atoms_p <- moldenATOMS
    
    -- parse the "[GTO]" block
    skipSpace
    basfuns_p <- moldenGTO
    
    -- parse the "[MO]" block
    skipSpace
    mmos_p <- moldenMO
    
    -- return the results
    return Molden { atoms = atoms_p
                  , basfuns = basfuns_p
                  , mos = mmos_p
                  }
    
moldenATOMS :: Parser (Units, [MoldenCoord])
moldenATOMS = do
    -- it starts with the "[Atoms]" block
    _ <- string $ T.pack "[Atoms]"
    
    -- parse the units. Maybe AU/(AU) or Angs/(Angs)
    skipSpace
    units_p <- do
        units_raw <- manyTill anyChar endOfLine
        if (T.isInfixOf (T.toLower . T.pack $ units_raw) (T.pack "(au)") == True)
           then return Bohr
           else return Angstrom
    
    -- the geometry line by line
    skipSpace
    geom <- many1 moldenCoordLineParser
    
    -- return the results
    return (units_p, geom)
    
    where
        -- parse a coordinate line in a [Atoms] block
        moldenCoordLineParser :: Parser MoldenCoord
        moldenCoordLineParser = do
            -- parse the atomic symbol
            skipSpace
            element_p <- manyTill anyChar (char ' ')
            
            -- parse the number of the atom in the molecule
            skipSpace
            natom_p <- decimal
            
            -- parse the atomic number (aka number of protons)
            skipSpace
            nproton_p <- decimal
            
            -- the x coordinate
            skipSpace
            x_p <- double
            
            -- the y coordinate
            skipSpace
            y_p <- double
            
            -- the z coordinate
            skipSpace
            z_p <- double
            
            --skipSpace
            -- return the results
            return MoldenCoord { moco_element = element_p
                               , moco_natom = natom_p
                               , moco_nproton = nproton_p
                               , moco_coord = (x_p, y_p, z_p)
                               }

moldenGTO :: Parser [[BasFun]]
moldenGTO = do
    -- it starts with the "[GTO]" block
    _ <- string $ T.pack "[GTO]"
    
    -- parse the basis functions of an atom
    skipSpace
    basfuns_p <- many1 moldenGTOAtomParser
    
    -- parse the informations about coordinate system
    skipSpace
    _ <- maybeOption (string $ T.pack "[5D]")
    skipSpace
    _ <- maybeOption (string $ T.pack "[7F]")
    skipSpace
    _ <- maybeOption (string $ T.pack "[9G]")
    skipSpace
    
    -- return the results (basis functions atomwise)
    return basfuns_p
    
    where
        {- parse the basis functions for a complete atom, for example
            2  0
          s      3  0
                13.0100000000         0.0334987264
                 1.9620000000         0.2348008012
                 0.4446000000         0.8136829579
          s      1  0
                 0.1220000000         1.0000000000
          p      1  0
                 0.7270000000         1.0000000000
        -}
        moldenGTOAtomParser :: Parser [BasFun]
        moldenGTOAtomParser = do
            -- parse the atom number, on which the basis functions are centred
            skipSpace
            _ <- (decimal :: Parser Int)
            
            -- the zero at the end
            skipSpace
            _ <- char '0'
            endOfLine
            
            -- parse the basis functions of the given atom
            basfunsAtom_p <- many1 moldenGTOAtomBFParser    
            
            --return the result
            return basfunsAtom_p
        
        
        {- parse a single basis function for a given atom 
          s      3  0
                13.0100000000         0.0334987264
                 1.9620000000         0.2348008012
                 0.4446000000         0.8136829579
        -}
        moldenGTOAtomBFParser :: Parser BasFun
        moldenGTOAtomBFParser = do
            -- parse the angular momentum of this basis function
            skipSpace
            angular_p <- do
                angular_char <- anyChar
                return $ orb2AngMom angular_char
            
            -- number of PGTOs to experct
            skipSpace
            npgto_p <- decimal
            
            -- parse the zero or other number till the endOfLine
            skipSpace
            _ <- manyTill anyChar endOfLine
            
            -- parse the PGTOs line by line and make a CGTO of them
            cgto_p <- count npgto_p pgto_and_ContrCoeff
            
            return $ BasFun { basfun_angular = angular_p
                            , basfun_radial = cgto_p
                            }
        
        -- parse a single pgto and its contraction coefficient
        pgto_and_ContrCoeff :: Parser (PGTO, ContrCoeff)
        pgto_and_ContrCoeff = do
            -- parse the exponent
            skipSpace
            pgto_p <- double
            
            -- parse the contraction coefficient
            skipSpace
            contrcoeff_p <- double
            endOfLine
            
            -- return the result
            return $ (pgto_p, contrcoeff_p)

moldenMO :: Parser [MMO]
moldenMO = do
    -- it starts with the "[Atoms]" block
    skipSpace
    _ <- string $ T.pack "[MO]"
    
    -- parse the MOs
    skipSpace
    mmos <- many1 singleMMOParser
    
    -- return the results
    return mmos
    
    where
        singleMMOParser :: Parser MMO
        singleMMOParser = do
            -- parse the "Sym" statement
            skipSpace
            _ <- string $ T.pack "Sym="
            skipSpace
            sym_p <- manyTill anyChar $ char ' '
            
            -- parse the "Ene" statement
            skipSpace
            _ <- string $ T.pack "Ene="
            skipSpace
            ene_p <- double
            
            -- parse the "Spin" statement
            skipSpace
            _ <- string $ T.pack "Spin="
            skipSpace
            spin_p <- do
                spin_raw <- manyTill anyChar $ (char ' ') <|> (char '\n')
                if (spin_raw == "Alpha")
                   then return Alpha
                   else return Beta
            
            -- parse the "Occup" statement
            skipSpace
            _ <- string $ T.pack "Occup="
            skipSpace
            occup_p <- double
            
            -- parse the MO coefficients in the AO expansion
            skipSpace
            coeffs_p <- do
                coeffs_raw <- many1 coeff_lineParser
                return $ BLAS.fromList coeffs_raw
            
            -- return the results
            return MMO { sym = sym_p
                       , energy = ene_p
                       , spin = spin_p
                       , occup = occup_p
                       , coeffs = coeffs_p 
                       }


        -- parse a single line of the coefficients in the [MO] block for a single MO
        coeff_lineParser :: Parser Double
        coeff_lineParser = do
            -- not interested in the first number (label of basis function)
            skipSpace            
            _ <- (decimal :: Parser Int)
            
            -- but parse the mo coefficient
            skipSpace
            single_mo_coeff_p <- double
            
            -- return the result
            return single_mo_coeff_p


{- ####################################### -}
{- Functions for working with molden files -}
{- ####################################### -}

-- get the overall number of electrons in the molecule
n_e :: Molden -> Double
n_e molden = sum . map occup . mos $ molden

-- get the overall number of basis functions in a molden file
n_BF :: Molden -> Int
n_BF molden = length . concat . basfuns$ molden

-- atomwise list of angular momentums for the basis functions
angMom_BF:: Molden -> [[Int]]
angMom_BF molden = map (map basfun_angular) $ basfuns molden

-- give in a basis function and you will get atomic orbitals from it in molden ordering
bf2aos_sph :: BasFun -> [AO]
bf2aos_sph bf
    -- this is the molden spherical ordering for s
    | bf_angmom == 0 = [ AO { ao_n = 0
                            , ao_l = bf_angmom
                            , ao_m = i
                            , ao_r = cgto
                            } | i <- [0]
                       ]
    -- this is the molden spherical ordering for p
    | bf_angmom == 1 = [ AO { ao_n = 0
                            , ao_l = bf_angmom
                            , ao_m = i
                            , ao_r = cgto
                            } | i <- [(-1), 0 , 1]
                       ]
    -- this is the molden spherical ordering for
    | bf_angmom == 2 = [ AO { ao_n = 0
                            , ao_l = bf_angmom
                            , ao_m = i
                            , ao_r = cgto
                            } | i <- [0, 1, (-1), 2, (-2)]
                       ]
    -- this is the molden spherical ordering for f
    | bf_angmom == 3 = [ AO { ao_n = 0
                            , ao_l = bf_angmom
                            , ao_m = i
                            , ao_r = cgto
                            } | i <- [0, 1, (-1), 2, (-2), 3, (-3)]
                       ]
    -- this is the molden spherical ordering for f
    | bf_angmom == 4 = [ AO { ao_n = 0
                            , ao_l = bf_angmom
                            , ao_m = i
                            , ao_r = cgto
                            } | i <- [0, 1, (-1), 2, (-2), 3, (-3), 4, (-4)]
                       ]
    | otherwise = error "Molden does not support functions higher than g"
    where bf_angmom = basfun_angular bf
          cgto = basfun_radial bf

-- from the basis functions expand in the AOs
-- ordering of the list
--   atom
--     Hauptquantenzahl, Nebequantenzahl (e.g. 1s / 2s / 3d / 4f / 5f)
--       Magnetquantenzahl
basis2AtomicOrbitals :: Molden -> [[[AO]]]
basis2AtomicOrbitals molden = mol_ao
    where
        -- the basis functions atomwise
        mol_bf = basfuns molden
        
        -- get angular momentums of the basis functions atomwise
        mol_ao = map (map bf2aos_sph) mol_bf

-- the order of the MOs in the basis of the AOs (giving angular momentum here)
moOrderingInAngMom :: Molden -> [Int]
moOrderingInAngMom molden = concat . concat . map (map (map ao_l)) . basis2AtomicOrbitals $ molden

-- the order of the MOs in the basis of the AOs (giving in corresponding atom)
moOrderingInAtoms :: Molden -> [Int]
moOrderingInAtoms molden =  fillLLwithPos . map concat . basis2AtomicOrbitals $ molden

-- the order of the MOs in the basis of the AOs (given in corresponding basis functions)
moOrderingInBF :: Molden -> [Int]
moOrderingInBF molden = fillLLwithPos . concat . map (map (map ao_l)) . basis2AtomicOrbitals $ molden

-- get MO coefficients (sum over AOs with same n and l but different m)
-- for a specific basis function in a molecular orbital
getMOcoeffsForBFinOrb :: Int -> Int -> Int -> Molden -> Double
getMOcoeffsForBFinOrb mo_number bf_number angMomOfBF molden
    | isInfixOf listIndexForBF listIndexForAngMom == True = sum coeffs_for_BF
    | otherwise = error "the desired angular momentum does not belong to the specified basis function"
    where
        mocoeffs = coeffs $ (mos molden) !! mo_number
        bfIndexList = moOrderingInBF molden
        angMomIndexList = moOrderingInAngMom molden
        listIndexForBF = findIndices (== bf_number) bfIndexList
        listIndexForAngMom = findIndices (== angMomOfBF) angMomIndexList
        coeffs_for_BF = [mocoeffs BLAS.! i | i <- listIndexForBF]

-- take a nested list and fills the inner entries with its position in the outer index, than flatten
-- [[1,2,3,4], [5,6,7]] -> [0,0,0,0,1,1,1]
fillLLwithPos :: [[a]] -> [Int]
fillLLwithPos a = concat [Prelude.take (length $ a !! i) $ repeat i | i <- [0 .. ((length a) - 1)]]

-- Make a parser optional, return Nothing if there is no match
maybeOption :: Parser a -> Parser (Maybe a)
maybeOption p = option Nothing (Just <$> p)
