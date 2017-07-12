module Data.Chemistry.XYZ
( XYZ(..)
, xyzParser
, getNAtoms
, getElements
, getCoords
, printXYZ
, getCoord2Mat
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
getCoord2Mat :: XYZ -> BLAS.Matrix Double
getCoord2Mat a = (nAtoms a BLAS.><3) $ concat . map (\(x, y, z) -> [x, y, z]) $ getCoords a
