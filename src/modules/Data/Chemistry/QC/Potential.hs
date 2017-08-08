module Data.Chemistry.QC.Potential
( cerjanMiller
, cerjanMiller'
, muellerBrown
, muellerBrown'
, vec_sin2D
, vec_sin2D'
, makeTestSurface
, print2D2Gnuplot
) where
import qualified Numeric.LinearAlgebra as BLAS
import System.IO
import Text.Printf

cerjanMiller :: Double -> Double -> Double -> BLAS.Vector Double -> Double
cerjanMiller a b c xy = (a - b * y**2) * x**2 * exp(-x**2) + c / 2 * y**2
    where
        x = xy BLAS.! 0
        y = xy BLAS.! 1

cerjanMiller' :: Double -> Double -> Double -> BLAS.Vector Double -> BLAS.Vector Double
cerjanMiller' a b c xy = 
    BLAS.vector [ 2 * (b * y**2 - a) * x**3 * exp(-x**2) - 2 * (b * y**2 - a) * x * exp(-x**2)
                , -2 * b * x**2 * y *exp(-x**2) + c * y ]
    where
        x = xy BLAS.! 0
        y = xy BLAS.! 1

muellerBrown :: ( BLAS.Vector Double
                , BLAS.Vector Double
                , BLAS.Vector Double
                , BLAS.Vector Double
                , BLAS.Vector Double
                , BLAS.Vector Double
                ) -> BLAS.Vector Double -> Double
muellerBrown (lAK, ak, bk, ck, x0k, y0k) xy = sum [ (lAK BLAS.! ind) * exp(
                                                                           (ak BLAS.! ind) * (x - (x0k BLAS.! ind))**2 +
                                                                           (bk BLAS.! ind) * (x - (x0k BLAS.! ind)) * (y - (y0k BLAS.! ind)) +
                                                                           (ck BLAS.! ind) * (y - (y0k BLAS.! ind))**2
                                                                          )
                                                  | ind <- [0 .. 3]
                                                  ]
    where
        x = xy BLAS.! 0
        y = xy BLAS.! 1


muellerBrown' :: ( BLAS.Vector Double
                 , BLAS.Vector Double
                 , BLAS.Vector Double
                 , BLAS.Vector Double
                 , BLAS.Vector Double
                 , BLAS.Vector Double
                 ) -> BLAS.Vector Double -> BLAS.Vector Double
muellerBrown' (lAK, ak, bk, ck, x0k, y0k) xy = BLAS.vector [x', y']
    where
        x' = 0.5 * ((muellerBrown (lAK, ak, bk, ck, x0k, y0k) xp_y) - (muellerBrown (lAK, ak, bk, ck, x0k, y0k) xy)) +
             0.5 * ((muellerBrown (lAK, ak, bk, ck, x0k, y0k) xy) - (muellerBrown (lAK, ak, bk, ck, x0k, y0k) xm_y))
        y' = 0.5 * ((muellerBrown (lAK, ak, bk, ck, x0k, y0k) x_yp) - (muellerBrown (lAK, ak, bk, ck, x0k, y0k) xy)) +
             0.5 * ((muellerBrown (lAK, ak, bk, ck, x0k, y0k) xy) - (muellerBrown (lAK, ak, bk, ck, x0k, y0k) x_ym))
        xp_y = xy + (BLAS.vector [0.05, 0.0])
        xm_y = xy + (BLAS.vector [-0.05, 0.0])
        x_yp = xy + (BLAS.vector [0.0, 0.05])
        x_ym = xy + (BLAS.vector [0.0, -0.05])

-- simple 2D sin surface
vec_sin2D :: BLAS.Vector Double -> Double
vec_sin2D xy = sin x + sin y
    where
        x = xy BLAS.! 0
        y = xy BLAS.! 1

vec_sin2D' :: BLAS.Vector Double -> BLAS.Vector Double
vec_sin2D' xy = BLAS.vector [cos x, cos y]
    where
        x = xy BLAS.! 0
        y = xy BLAS.! 1        

-- create a good miller brown surface
-- creates a list of tuples in the form of (vector coordinates, energy)
makeTestSurface :: [(BLAS.Vector Double, Double)]
makeTestSurface = zip xy $ map (muellerBrown (lAK, ak, bk, ck, x0k, y0k)) xy 
    where
        lAK = BLAS.vector [-200, -100, -170, 15]
        ak = BLAS.vector [-1, -1, -6.5, 0.7]
        bk = BLAS.vector [0, 0, 11, 0.6]
        ck = BLAS.vector [-10, -10, -6.5, 0.7]
        x0k = BLAS.vector [1, 0, -0.5, -1]
        y0k = BLAS.vector [0, 0.5, 1.5, 1]
        xy = map BLAS.fromList [[x, y] | x <- [-3.0, -2.95 .. 3.0], y <- [-3.0, -2.95 .. 3.0]]

-- print coordinate / energy pairs to gnuplot compatible format
print2D2Gnuplot :: Handle -> [(BLAS.Vector Double, Double)] -> IO ()
print2D2Gnuplot handle inCoordsValues = do
    mapM_ (print2DVal2Gnuplot handle) inCoordsValues

print2DVal2Gnuplot :: Handle -> (BLAS.Vector Double, Double) -> IO ()
print2DVal2Gnuplot handle (inCoord,value) = do
    hPrintf handle "%16.8F    " x
    hPrintf handle "%16.8F    " y
    hPrintf handle "%16.8F\n" value
    where
        x = inCoord BLAS.! 0
        y = inCoord BLAS.! 1
