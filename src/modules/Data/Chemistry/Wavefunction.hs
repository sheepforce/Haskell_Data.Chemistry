module Data.Chemistry.Wavefunction
( Spin(..)
, PGTO
, ContrCoeff
, CGTO
, BasFun(..)
, AO(..)
, MO
) where
import qualified Numeric.LinearAlgebra as BLAS
--import Data.Maybe

{- ################ -}
{- basic Data Types -}
{- ################ -}

-- spin of an electron may be alpha or beta
data Spin = Alpha | Beta deriving Show

-- a primitive GTO (aka a single Gaussian) is just an exponent
type PGTO = Double

-- a contraction coefficients is just a number
type ContrCoeff = Double

-- a contracted GTO is a weighted sum of Gaussians (coefficients and PGTOs)
type CGTO = [(PGTO,  ContrCoeff)]

-- a basis function is described by the nebenquantenzahl l (angular momentum)
-- and a radial function (CGTO)
data BasFun = BasFun { basfun_angular :: Int
                     , basfun_radial  :: CGTO
                     } deriving Show

-- an atomic orbital is characterized by
--   n = Hauptquantenzahl
--   l = Nebenquantenzahl
--   m = Magnetquantenzahl
--   m_s = Spinquantenzahl
--   the radial part from a CGTO
data AO = AO { ao_n :: Int
             , ao_l :: Int
             , ao_m :: Int
             , ao_r :: CGTO
             } deriving Show

-- a molecular orbital can be described by MO coefficients (LCAO coefficients)
type MO = BLAS.Vector Double


{- ############ -}
{- calculations -}
{- ############ -}
-- give a list of weighted MOs and get back a new MO built of these and renormalized
mixMOs :: [(Double, MO)] -> MO
mixMOs lcl = newWMO
  where
    norm = sum . map fst $ lcl
    normCoeffs = map (/ norm). map fst $ lcl :: [Double]
    mos = map snd lcl :: [MO]
    rlcl = zip normCoeffs mos :: [(Double, MO)]
    nBF = length . BLAS.toList . head $ mos
    zeroMO = BLAS.fromList $ replicate nBF 0.0 :: MO
    newWMO = foldr (\(c, m) acc -> BLAS.fromList [c] * m + acc) zeroMO rlcl
