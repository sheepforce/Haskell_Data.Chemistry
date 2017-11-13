module Data.Chemistry.Molden
( moldenParser
) where
import qualified Numeric.LinearAlgebra as BLAS
import qualified Data.Text as T
import Control.Applicative
import Data.Attoparsec.Text.Lazy
--import Data.Chemistry.XYZ
import Data.Chemistry.Wavefunction
import Data.Chemistry.BasisSet

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
    skipSpace
    
    -- parse the "[Atoms]" block, containing the geometry informations
    atoms_p <- do
        -- the beginning of the block is "[Atoms"]
        _ <- string $ T.pack "[Atoms]"
        skipSpace
        -- the units used
        units_p <- do
            units_raw <- (string $ T.pack "AU") <|> (string $ T.pack "Angs")
            if (units_raw == T.pack "AU")
               then return Bohr
               else return Angstrom
        -- the geometry line by line
        geom <- many1 moldenCoordLineParser
        return (units_p, geom)
    skipSpace
    
    -- parse the "[GTO]" block
    _ <- string $ T.pack "[GTO]"
    skipSpace
    -- parse (number of atoms) basises
    basfuns_p <- count (length . snd $ atoms_p) moldenGTOAtomParser
    skipSpace
    
    -- parse the "[MO]" block
    _ <- string $ T.pack "[MO]"
    skipSpace
    mmos_p <- many1 moldenMMOParser 
    
    return Molden { atoms = atoms_p
                  , basfuns = basfuns_p
                  , mos = mmos_p
                  }
    
    


-- parse a coordinate line in a [Atoms] block
moldenCoordLineParser :: Parser MoldenCoord
moldenCoordLineParser = do
    skipSpace
    element_p <- manyTill anyChar (char ' ')
    skipSpace
    natom_p <- decimal
    skipSpace
    nproton_p <- decimal
    skipSpace
    x_p <- double
    skipSpace
    y_p <- double
    skipSpace
    z_p <- double
    skipSpace
    return MoldenCoord { moco_element = element_p
                       , moco_natom = natom_p
                       , moco_nproton = nproton_p
                       , moco_coord = (x_p, y_p, z_p)
                       }

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
    basfuns_p <- many1 moldenGTOAtomBFParser    
    
    return basfuns_p

{- parse a single basis function for a given atom 
  s      3  0
        13.0100000000         0.0334987264
         1.9620000000         0.2348008012
         0.4446000000         0.8136829579
-}
moldenGTOAtomBFParser :: Parser BasFun
moldenGTOAtomBFParser = do
    skipSpace
    -- parse the angular momentum of this basis function
    angular_p <- do
        angular_char <- anyChar
        return $ orb2AngMom angular_char
    
    -- number of PGTOs to experct
    skipSpace
    npgto_p <- decimal
    
    -- parse the zero
    skipSpace
    _ <- char '0'
    endOfLine
    
    -- parse the PGTOs line by line and make a CGTO of them
    cgto_p <- count npgto_p pgto_and_ContrCoeff
    
    return $ BasFun { basfun_angular = angular_p
                    , basfun_radial = cgto_p
                    }
    
    where
        -- parse a single pgto and its contraction coefficient
        pgto_and_ContrCoeff :: Parser (PGTO, ContrCoeff)
        pgto_and_ContrCoeff = do
            skipSpace
            pgto_p <- double
            skipSpace
            contrcoeff_p <- double
            endOfLine
            return $ (pgto_p, contrcoeff_p)
    
moldenMMOParser :: Parser MMO
moldenMMOParser = do
    -- parse the "Sym" statement
    _ <- string $ T.pack "Sym="
    skipSpace
    sym_p <- manyTill anyChar $ char ' '
    skipSpace
    
    -- parse the "Ene" statement
    _ <- string $ T.pack "Ene="
    skipSpace
    ene_p <- double
    skipSpace
    
    -- parse the "Spin" statement
    _ <- string $ T.pack "Spin="
    skipSpace
    spin_p <- do
        spin_raw <- manyTill anyChar $ (char ' ') <|> (char '\n')
        if (spin_raw == "Alpha")
           then return Alpha
           else return Beta
    skipSpace
    
    -- parse the "Occup" statement
    _ <- string $ T.pack "Occup="
    skipSpace
    occup_p <- double
    skipSpace
    
    -- parse the MO coefficients in the AO expansion
    coeffs_p <- do
        coeffs_raw <- many1 coeff_lineParser
        return $ BLAS.fromList coeffs_raw
    
    skipSpace
    
    return MMO { sym = sym_p
               , energy = ene_p
               , spin = spin_p
               , occup = occup_p
               , coeffs = coeffs_p 
               }
    where
        -- parse a single line of the coefficients in the [MO] block for a single MO
        coeff_lineParser :: Parser Double
        coeff_lineParser = do
            skipSpace
            -- not interested in the number of bs used
            _ <- double
            skipSpace
            single_mo_coeff_p <- double
            return single_mo_coeff_p
