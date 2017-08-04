module Data.Chemistry.QC
( QcInMolecule(..)
) where
import qualified Numeric.LinearAlgebra as BLAS
--import qualified Data.Chemistry.XYZ as XYZ

data QcInMolecule = QcInMolecule { geom_cart :: BLAS.Vector Double
                                 , multiplicity :: Int
                                 , charge :: Int
                                 }
