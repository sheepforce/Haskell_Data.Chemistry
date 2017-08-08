module Data.Chemistry.QC
( QcInMolecule(..)
, interpolateGeom
, vectorProjection
, steepestDescent
, bfgs
, bfgs2
) where
import qualified Numeric.LinearAlgebra as BLAS
import qualified Numeric.GSL as GSL
--import qualified Data.Chemistry.XYZ as XYZ


{-
EXAMPLES:
    1) -- optimize to the maximum on the 2D sin surface from start geometry at [-0.5, -0.5] and end geometry at [-0.93, 3.47]
       let start = BLAS.vector [-0.5, -0.5]
       let end = BLAS.vector [-0.93, 3.47]
       let reactionVector = end - start
       steepestDescent 1.0e-6 1000 0.05 1.0e-5 (negate . vec_sin2D) (negate. (vectorProjection reactionVector) . vec_sin2D') (BLAS.vector [-0.59, 0.41])
       (correct result at [-0.706796744008091, 1.4883327295630742])
    2) -- optimize to the maximum on the path between Müller-Brown surface with start at [0.6211317901475866,  2.823674628779631e-2]
       -- and end at [-0.5591208056179865, 1.4407836667987057]
       let start = BLAS.vector [0.6211317901475866,  2.823674628779631e-2]
       let end = BLAS.vector [-0.5591208056179865, 1.4407836667987057]
       let reaction_vec = end - start
       steepestDescent 1.0e-6 1000 0.02 1.0e-5 (negate . muellerBrown (lAK, ak, bk, ck, x0k, y0k)) (negate . (vectorProjection reaction_vec) . muellerBrown' (lAK, ak, bk, ck, x0k, y0k)) end
       (correct result at [-0.19720416827758844, 1.0076355141142914])
-}


data QcInMolecule = QcInMolecule { geom_cart :: BLAS.Vector Double
                                 , multiplicity :: Int
                                 , charge :: Int
                                 }

-- take two geometries and interpolate on a linear path between them
interpolateGeom :: (BLAS.Vector Double, BLAS.Vector Double) -> Int -> [BLAS.Vector Double]
interpolateGeom (geom1, geom2) nImgs = [(BLAS.vector [p])  * geom2 + (1 - (BLAS.vector [p])) * geom1 | p <- (BLAS.toList $ BLAS.linspace nImgs (0, 1 :: Double))]

-- project one vector on an other (a on b)
vectorProjection :: BLAS.Vector Double -> BLAS.Vector Double -> BLAS.Vector Double
vectorProjection b_vec a_vec = (BLAS.vector [a1]) * bUnit_vec
    where
        a1 = a_vec BLAS.<.> bUnit_vec
        bUnit_vec = b_vec / (BLAS.fromList [(sqrt . BLAS.sumElements $ (b_vec * b_vec))])

-- perform a steepest descent optimization with gradients as implemented in GSL/hmatrix-gsl
-- to make it an ascending algorithm, negate the energy and gradient and project the gradient 
-- on the desired vector
steepestDescent ::
    Double ->                                       -- desired precision of the solution (gradient test)
    Int ->                                          -- maximum number of iterations
    Double ->                                       -- size of first trial step
    Double ->                                       -- tolerance - see GSL manual, recommended 0.1, sets step size
    (BLAS.Vector Double -> Double) ->               -- function to optimize
    (BLAS.Vector Double -> BLAS.Vector Double) ->   -- gradient of the function
    BLAS.Vector Double ->                           -- starting point
    (BLAS.Vector Double, BLAS.Matrix Double)        -- solution vector and optimization path
steepestDescent gradConv maxIter trialStepSize tolerance energy gradient initGeom = 
    GSL.minimizeVD GSL.SteepestDescent gradConv maxIter trialStepSize tolerance energy gradient initGeom

-- BFGS optimization algorithm, also works for the projected maximum when negating energy and gradient and
-- projecting the gradient
bfgs ::
    Double ->                                       -- desired precision of the solution (gradient test)
    Int ->                                          -- maximum number of iterations
    Double ->                                       -- size of first trial step
    Double ->                                       -- tolerance - see GSL manual, recommended 0.1, sets step size
    (BLAS.Vector Double -> Double) ->               -- function to optimize
    (BLAS.Vector Double -> BLAS.Vector Double) ->   -- gradient of the function
    BLAS.Vector Double ->                           -- starting point
    (BLAS.Vector Double, BLAS.Matrix Double)        -- solution vector and optimization path
bfgs gradConv maxIter trialStepSize tolerance energy gradient initGeom = 
    GSL.minimizeVD GSL.VectorBFGS gradConv maxIter trialStepSize tolerance energy gradient initGeom

-- BFGS optimization algorithm, also works for the projected maximum when negating energy and gradient and
-- projecting the gradient
-- RECOMMENDED ALGORITHM
bfgs2 ::
    Double ->                                       -- desired precision of the solution (gradient test)
    Int ->                                          -- maximum number of iterations
    Double ->                                       -- size of first trial step
    Double ->                                       -- tolerance - see GSL manual, recommended 0.1, sets step size
    (BLAS.Vector Double -> Double) ->               -- function to optimize
    (BLAS.Vector Double -> BLAS.Vector Double) ->   -- gradient of the function
    BLAS.Vector Double ->                           -- starting point
    (BLAS.Vector Double, BLAS.Matrix Double)        -- solution vector and optimization path
bfgs2 gradConv maxIter trialStepSize tolerance energy gradient initGeom = 
    GSL.minimizeVD GSL.VectorBFGS gradConv maxIter trialStepSize tolerance energy gradient initGeom




   
