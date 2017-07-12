module Data.Chemistry.XYZ
( XYZ(..)
, xyzParser
, getNAtoms
, getElements
, getCoords
, printXYZ
, coord2Mat
, mat2Coord
, interpolate
) where
import System.IO
import Data.Attoparsec.Text.Lazy
import Text.Printf
import qualified Numeric.LinearAlgebra as BLAS


{- ################# -}
{- define Data Types -}
{- ################# -}

data XYZ = XYZ { nAtoms :: Int
               , comment :: String
               , xyzcontent :: [(String,Double,Double,Double)]
               } deriving Show

{- ###### -}
{- Parser -}
{- ###### -}

xyzParser :: Parser XYZ
xyzParser = do
  _ <- many' $ char ' '
  nAtoms_parse <- decimal
  _ <- many' $ char ' '
  _ <- endOfLine
  comment_parse <- manyTill anyChar endOfLine
  -- coordinates <- many' xyzCoordLineParser <* (endOfLine <|> endOfInput)
  coordinates <- many' xyzCoordLineParser
  _ <- endOfInput
  return $ XYZ { nAtoms = nAtoms_parse
               , comment = comment_parse
               , xyzcontent = coordinates
               }
    where
      xyzCoordLineParser :: Parser (String,Double,Double,Double)
      xyzCoordLineParser = do
        _ <- many' $ char ' '        
        element <- manyTill anyChar (char ' ')
        _ <- many' $ char ' '
        x <- double
        _ <- many' $ char ' '
        y <- double
        _ <- many' $ char ' '
        z <- double
        _ <- many' $ char ' '
        _ <- many' endOfLine
        return $ (element,x,y,z)


{- ################################ -}
{- Functions to work with XYZ files -}
{- ################################ -}

getNAtoms :: XYZ -> Int
getNAtoms a = nAtoms a

getElement :: (String,Double,Double,Double) -> String
getElement (e, _, _, _) = e

getElements :: XYZ -> [String]
getElements a = map getElement $ xyzcontent a

getCoord :: (String,Double,Double,Double) -> (Double,Double,Double)
getCoord (_, x, y, z) = (x, y, z)

getCoords :: XYZ -> [(Double,Double,Double)]
getCoords a = map getCoord $ xyzcontent a

-- print XYZ formated to a handle
-- using mapM_ over a list enables printing a trajectory
printXYZ :: Handle -> XYZ ->  IO ()
printXYZ xyzhandle xyz = do
    hPrint xyzhandle $ nAtoms xyz
    hPrint xyzhandle $ comment xyz
    mapM_ (\(e, x, y, z) -> do
        hPrintf xyzhandle "%-4s" e
        hPrintf xyzhandle "%+16.8f    " x
        hPrintf xyzhandle "%+16.8f    " y
        hPrintf xyzhandle "%+16.8f\n" z) (xyzcontent xyz)

-- generate HMatrix Matrix of coordinates
coord2Mat :: XYZ -> BLAS.Matrix Double
coord2Mat a = (nAtoms a BLAS.>< 3) $ concat . map (\(x, y, z) -> [x, y, z]) $ getCoords a

-- fill XYZ with coordinates from a matrix (replacing the old ones)
mat2Coord :: XYZ -> BLAS.Matrix Double -> XYZ
mat2Coord a mat = XYZ { nAtoms = (nAtoms a)
                      , comment = (comment a)
                      , xyzcontent = content
                      }
    where
        content = zipWith (\e (x, y, z) -> (e, x, y, z)) (getElements a) coords
        coords = (listList2TupleList . mat2ListList) mat
        
        mat2ListList :: BLAS.Matrix Double -> [[Double]]
        mat2ListList mat_in = [BLAS.toList (mat_in BLAS.! i) | i <- [0 .. (BLAS.rows mat_in - 1)]]
        
        listList2TupleList :: [[Double]] -> [(Double, Double, Double)]
        listList2TupleList listlist = map (\[x,y,z] -> (x,y,z)) listlist

-- interpolate two geometries in cartesian coordinates
interpolate :: Int -> XYZ -> XYZ -> [XYZ]
interpolate nImages mol1 mol2 = map (mat2Coord mol1) interpol_mats
    where
        interpol_mats = [mol1_coordMat * ((1 BLAS.>< 1) [(1.0 - p)]) + mol2_coordMat * ((1 BLAS.>< 1) [p]) | p <- (BLAS.toList $ BLAS.linspace nImages (0, 1 :: Double))]
        mol1_coordMat = coord2Mat mol1
        mol2_coordMat = coord2Mat mol2
