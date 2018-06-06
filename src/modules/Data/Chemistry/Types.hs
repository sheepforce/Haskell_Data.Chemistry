{-# LANGUAGE TemplateHaskell #-}
module Data.Chemistry.Types
( Spin(..)
, PGTO
, ContrCoeff
, CGTO
, BasFun(..)
, basfun_angular
, basfun_radial
, AO(..)
, ao_n
, ao_l
, ao_m
, ao_r
, MO
, ConGauss(..)
, Basis(..)
, XYZ(..)
, xyz_nAtoms
, xyz_comment
, xyz_xyzcontent
, MoldenUnits(..)
, MoldenCoord(..)
, moldenCoord_element
, moldenCoord_natom
, moldenCoord_nproton
, moldenCoord_coord
, MoldenMO(..)
, moldenMO_sym
, moldenMO_energy
, moldenMO_spin
, moldenMO_occup
, moldenMO_coeffs
, Molden(..)
, molden_atoms
, molden_basfuns
, molden_mos
) where
import           Lens.Micro.Platform
import qualified Numeric.LinearAlgebra as BLAS


--------------------------------------------------------------------------------
-- general abstract handling of Wavefunctions
--------------------------------------------------------------------------------
-- | Spin of an electron may be alpha or beta
data Spin = Alpha | Beta deriving Show

-- | A primitive GTO (aka a single Gaussian) is just an exponent
type PGTO = Double

-- | A contraction coefficients is just a number
type ContrCoeff = Double

-- | A contracted GTO is a weighted sum of Gaussians (coefficients and PGTOs)
type CGTO = [(PGTO,  ContrCoeff)]

-- | A basis function is described by the nebenquantenzahl l (angular momentum)
-- | and a radial function (CGTO)
data BasFun = BasFun
  { _basfun_angular :: Int
  , _basfun_radial  :: CGTO
  } deriving Show
makeLenses ''BasFun

-- | An atomic orbital is characterized by
-- |   n = Hauptquantenzahl
-- |   l = Nebenquantenzahl
-- |   m = Magnetquantenzahl
-- |   m_s = Spinquantenzahl
-- |   the radial part from a CGTO
data AO = AO
  { _ao_n :: Int
  , _ao_l :: Int
  , _ao_m :: Int
  , _ao_r :: CGTO
  } deriving Show
makeLenses ''AO

-- | A molecular orbital can be described by MO coefficients (LCAO coefficients)
type MO = BLAS.Vector Double

-- | A single contracted gaussian made of primitve gaussians
-- | (Double, Double) contains exponents and coefficients
data ConGauss = ConGauss
  { angMom          :: Int
  , expoCoeff_pairs :: [(Double, Double)]
  } deriving Show

-- | Define data representing a finite gaussian basis
data Basis = Basis
  { element :: String
  , basFuns :: [ConGauss]
  } deriving Show


--------------------------------------------------------------------------------
-- chemical data formats
--------------------------------------------------------------------------------
{- === -}
{- XYZ -}
{- === -}
-- | contents of XYZ files (no connectivities)
data XYZ = XYZ
  { _xyz_nAtoms     :: Int
  , _xyz_comment    :: String
  , _xyz_xyzcontent :: [(String,Double,Double,Double)]
  } deriving Show
makeLenses ''XYZ

{- ====== -}
{- Molden -}
{- ====== -}
-- | units for geometry may either be Angstrom or Bohr
data MoldenUnits = Angstrom | Bohr deriving (Show, Eq)

-- | A data type for the geometry informations in a molden [Atoms] block
data MoldenCoord = MoldenCoord
  { _moldenCoord_element :: String
  , _moldenCoord_natom   :: Int
  , _moldenCoord_nproton :: Int
  , _moldenCoord_coord   :: (Double, Double, Double)
  } deriving Show
makeLenses ''MoldenCoord

-- | Extended data type for MOs in a Molden-File
data MoldenMO = MoldenMO
  { _moldenMO_sym    :: String
  , _moldenMO_energy :: Double
  , _moldenMO_spin   :: Spin
  , _moldenMO_occup  :: Double
  , _moldenMO_coeffs :: MO
  } deriving Show
makeLenses ''MoldenMO

-- | Contents of a Molden file
data Molden = Molden
  { _molden_atoms   :: (MoldenUnits, [MoldenCoord]) -- atoms contains the geometry
  , _molden_basfuns :: [[BasFun]]                   -- the basis functions, outer for the atom, inner for gto
  , _molden_mos     :: [MoldenMO]
  } deriving Show
makeLenses ''Molden
